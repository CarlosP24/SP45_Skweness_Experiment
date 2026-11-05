# Step 1 of the devices-B/C fit: digitize the published SI figures (first estimate, until the lab files
# of devices B and C are available). Run from the SP45 root:
#     julia --project=plots refit_BC/digitize.jl [path/to/Overleaf_repo]
# Needs qpdf. Writes data_exp/digitized/*.dat and refit_BC/figures/digitization_check.png.
#
# Output format (DelimitedFiles, '#' comment header):
#   *_dVdI_*.dat : matrix. Row 1 = B (T) of each column (first entry NaN); column 1 = I (nA) of each row;
#                  entries dV/dI (Ω) with the published R_s already subtracted. One row per image pixel
#                  row (the published maps are interpolated in I, so this oversamples the measurement);
#                  columns on the field grid of the measurement (from the R_0 polyline of panel c).
#   *_R0_line.dat: two columns, B (T) and the published zero-bias R_0 (Ω) of panel (c).
using DelimitedFiles, Printf, Statistics, CairoMakie
include(joinpath(pwd(), "plots", "kitchen", "pdf_digitizer.jl"))
using .PDFDigitizer

const OVERLEAF = isempty(ARGS) ?
    joinpath(homedir(), "Library/CloudStorage/Nextcloud-saco.csic.es-CARLOSPAYAHERRERO/PhD/Preprints/2026_07_SkewedIc/Overleaf_repo") :
    ARGS[1]
const OUT = joinpath("data_exp", "digitized")
const FIG = joinpath("refit_BC", "figures")
mkpath(OUT); mkpath(FIG)

const C0BLUE = (0.122, 0.467, 0.706)      # matplotlib C0, the R_0(B) line of panel (c)

# Tick values read off the rendered figures (the labels are outlined glyphs, positions are fitted)
const SPECS = Dict(
    "B" => (pdf = "SI_DEVICE_B.pdf", Rs = 5490.0, Bstep = 0.004, Bmax = 0.24,
            line = (xt = 0:0.05:0.25, yt = 0:200:1000),
            maps = [(tag = "a", W = 2034, cbW = 104, xt = 0:0.05:0.20, yt = 1e3 .* (-0.6:0.2:0.6), cbt = 0:250:2000)]),
    "C" => (pdf = "SI_DEVICE_C.pdf", Rs = 5350.0, Bstep = 0.003, Bmax = 0.48,
            line = (xt = 0:0.1:0.5, yt = 0:200:1000),
            maps = [(tag = "b", W = 2065, cbW = 89, xt = 0:0.1:0.4, yt = -75:25:100, cbt = 0:250:2000),
                    (tag = "a", W = 2043, cbW = 88, xt = 0:0.1:0.4, yt = 1e3 .* (-0.8:0.2:0.6), cbt = 0:250:2000)]),
)

"R_0(B) polyline of panel (c), in data units, and the residual of the tick fits"
function read_line(paths, spec)
    p = only(strokes_in(paths, (-Inf, Inf, -Inf, Inf), C0BLUE))
    pts = polyline(p)
    fr = frame_of(paths, pts[1][1] + 1, pts[1][2] + 1)
    fx, rx = fit_axis(ticks_below(paths, fr[3], fr[1:2]), collect(spec.xt))
    fy, ry = fit_axis(ticks_left(paths, fr[1], fr[3:4]), collect(spec.yt))
    return [fx(q[1]) for q in pts], [fy(q[2]) for q in pts], (rx, ry)
end

"""
Invert one dV/dI map. Returns B per pixel column, I (nA) per pixel row (top row first), the values,
and diagnostics (tick-fit residuals, colours off the colourbar, saturated fraction, colour level).
"""
function read_map(objs, paths, images, m)
    im = only([i for i in images if i.W == m.W])
    cb = only([i for i in images if i.W == m.cbW])
    xs, ys = pixel_centres(im)
    fr = frame_of(paths, xs[end ÷ 2], ys[end ÷ 2])
    fx, rx = fit_axis(ticks_below(paths, fr[3], fr[1:2]), collect(m.xt))
    fy, ry = fit_axis(ticks_left(paths, fr[1], fr[3:4]), collect(m.yt))
    e = image_extent(cb)
    fcb, rcb = fit_axis(ticks_right(paths, e[2], (e[3] - 0.2, e[4] + 0.2); tol = 0.1), collect(m.cbt))
    lut, levels = colorbar_lut(image_pixels(objs, cb.id), cb, fcb)
    vals, nmiss = invert_image(image_pixels(objs, im.id), lut, levels)
    top = levels[end][2]
    diag = (tick_res = (rx, ry, rcb), nlevels = length(levels), level = median(diff(last.(levels))),
            offlut = nmiss / length(vals), saturated = count(>=(top), vals) / length(vals), top = top)
    return fx.(xs), fy.(ys), vals, diag
end

"linear interpolation of the map columns at the measurement fields Bgrid"
function sample_columns(Bpix, vals, Bgrid)
    out = Matrix{Float64}(undef, size(vals, 1), length(Bgrid))
    for (k, b) in enumerate(Bgrid)
        j = clamp(searchsortedlast(Bpix, b), 1, length(Bpix) - 1)
        t = clamp((b - Bpix[j]) / (Bpix[j+1] - Bpix[j]), 0, 1)
        out[:, k] .= (1 - t) .* vals[:, j] .+ t .* vals[:, j+1]
    end
    return out
end

function write_map(path, Bgrid, Is, M, header)
    open(path, "w") do io
        for h in header
            println(io, "# ", h)
        end
        writedlm(io, permutedims(vcat(NaN, Bgrid)))
        # rows from the most negative current up
        for i in sortperm(Is)
            writedlm(io, permutedims(vcat(round(Is[i]; digits = 4), round.(M[i, :]; digits = 2))))
        end
    end
end

report = String[]
checks = Dict{String,Any}()
for dev in ("B", "C")
    sp = SPECS[dev]
    objs, paths, images = digitize_pdf(joinpath(OVERLEAF, sp.pdf))
    Bl, R0l, lres = read_line(paths, sp.line)
    Bgrid = collect(0:sp.Bstep:sp.Bmax)
    # the polyline is matplotlib-simplified (collinear points dropped): check it sits on the grid
    off = maximum(abs.(Bl ./ sp.Bstep .- round.(Bl ./ sp.Bstep)))
    off < 0.01 || error("device $dev: R_0 line not on a $(sp.Bstep) T grid (worst offset $off steps)")
    open(joinpath(OUT, "dev$(dev)_R0_line.dat"), "w") do io
        println(io, "# Device $dev, $(sp.pdf) panel (c): published zero-bias R_0(B), R_s = $(sp.Rs) Ω subtracted")
        println(io, "# columns: B (T), R0 (Ω). Vector polyline, $(length(Bl)) points (path-simplified on a $(sp.Bstep) T grid)")
        writedlm(io, hcat(round.(Bl; digits = 6), round.(R0l; digits = 3)))
    end
    push!(report, @sprintf("device %s: R_0 line %d points on the %.3f T grid, tick-fit residuals %.1e T, %.1e Ω",
        dev, length(Bl), sp.Bstep, lres...))
    for m in sp.maps
        Bpix, Ipix, vals, d = read_map(objs, paths, images, m)
        M = sample_columns(Bpix, vals, Bgrid)
        kind = m.tag == "b" || dev == "B" ? (dev == "B" ? "map" : "zoom") : "full"
        path = joinpath(OUT, "dev$(dev)_dVdI_$(kind).dat")
        write_map(path, Bgrid, Ipix, M, [
            "Device $dev, $(sp.pdf) panel ($(m.tag)): dV/dI (Ω) after the published R_s = $(sp.Rs) Ω subtraction",
            "digitized by refit_BC/digitize.jl (plots/kitchen/pdf_digitizer.jl): colour → value through the embedded colourbar",
            @sprintf("colour level %.2f Ω (256 levels over %d–%d Ω); values ≥ %.1f Ω are clipped (%.2f %% of pixels)",
                d.level, first(m.cbt), last(m.cbt), d.top, 100d.saturated),
            "row 1: B (T) of each column; column 1: I (nA) of each row (one row per image pixel row, oversampled)",
        ])
        push!(report, @sprintf("device %s panel (%s) → %s: %d × %d, B %.4f–%.4f T, I %.1f–%.1f nA, level %.2f Ω, off-LUT colours %.1f %% (all ±1 rounding), clipped %.2f %%, tick residuals %.1e T, %.1e nA, %.1e Ω",
            dev, m.tag, basename(path), length(Ipix), length(Bgrid), Bpix[1], Bpix[end], minimum(Ipix), maximum(Ipix),
            d.level, 100d.offlut, 100d.saturated, d.tick_res...))
        # check against the published R_0 line: the map row at I ≈ 0, interpolated at the line's fields
        i0 = argmin(abs.(Ipix))
        r0map = [M[i0, argmin(abs.(Bgrid .- b))] for b in Bl]
        dev_ = r0map .- R0l
        push!(report, @sprintf("   map at I = %.2f nA vs the R_0 line: mean %+.1f Ω, rms %.1f Ω, max |Δ| %.1f Ω",
            Ipix[i0], mean(dev_), sqrt(mean(dev_ .^ 2)), maximum(abs.(dev_))))
        checks["$dev$(m.tag)"] = (; Bgrid, Ipix, M, Bl, R0l, r0map, i0)
    end
end
open(joinpath(OUT, "digitization_log.txt"), "w") do io
    foreach(l -> println(io, l), report)
end
foreach(println, report)

# Validation figure: the digitized maps (same colour scale as the SI) and the I ≈ 0 cut vs the line
fig = Figure(size = (1500, 900), fontsize = 14)
for (k, key) in enumerate(("Ba", "Cb", "Ca"))
    c = checks[key]
    ax = Axis(fig[1, k]; title = "digitized $(key[1]) panel ($(key[2]))", xlabel = "B (T)", ylabel = "I (nA)")
    p = sortperm(c.Ipix)
    hm = heatmap!(ax, c.Bgrid, c.Ipix[p], permutedims(c.M[p, :]); colormap = :seismic, colorrange = (0, 2000))
    k == 3 && Colorbar(fig[1, 4], hm; label = "dV/dI (Ω)")
    ax2 = Axis(fig[2, k]; title = "R_0: published line (black) vs map at I ≈ 0 (orange)", xlabel = "B (T)", ylabel = "R_0 (Ω)")
    lines!(ax2, c.Bl, c.R0l; color = :black, linewidth = 2)
    scatterlines!(ax2, c.Bl, c.r0map; color = :orange, markersize = 4)
end
save(joinpath(FIG, "digitization_check.png"), fig)
println("saved ", joinpath(FIG, "digitization_check.png"))
