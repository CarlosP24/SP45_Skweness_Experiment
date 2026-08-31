# Hero slot — 771 × 240 mm. Paper Fig. 3, rebuilt wide as a 2 × 2 grid.
#
#   ┌──────────────────────┬──────────────────────┐
#   │ (a) LDOS, left lead   │ (c) dV/dI + Ic       │  ~60 % of right column
#   ├──────────────────────┤                      │
#   │ (b) LDOS, right lead  ├──────────────────────┤
#   │                      │ (d) R₀/R_N            │  ~40 %
#   └──────────────────────┴──────────────────────┘
#
# Re-style only: every number comes from data/ and data_exp/ through the paper's
# own `process_data`, which is included unmodified from plots/plotters/.
#
# Run from anywhere:  julia --project=plots poster/hero_fig3.jl

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)   # the paper's readers use paths relative to the repo root

using Parameters, JLD2, Quantica, FullShell, DelimitedFiles
using CairoMakie, Statistics, SavitzkyGolay, Roots, Bessels, Peaks

include(joinpath(ROOT, "poster", "style", "poster_theme.jl"))
include(joinpath(ROOT, "plots", "kitchen", "data_reader.jl"))
include(joinpath(ROOT, "plots", "kitchen", "Ic_reader.jl"))
include(joinpath(ROOT, "plots", "plotters", "theo_exp_core.jl"))  # process_data, ics — unmodified

set_theme!(poster_theme)

const BLIMS      = (0.0, 0.4)        # T — as in the paper, not one gauss wider
const LDOS_RANGE = (0, 2e-3)
const DVDI_RANGE = (5e2, 1.5e3)      # Ω
const CMAP_C     = :seismic          # the paper's own dV/dI map, kept as published

"LDOS panel data, loaded exactly as the paper's `panel_LDOS` does."
function load_ldos(name::String, side::String; basepath = "data/LDOS")
    res = load(joinpath(basepath, "fitted_$(name)_$(side).jld2"))["res"]
    @unpack LDOS, params = res
    @unpack Brng = params
    ωrng = real.(collect(params.ωrng))
    if side == "left"                       # mirror the left lead onto ω > 0
        LDOS = reverse(LDOS, dims = 2)
        ωrng = -reverse(ωrng)
    end
    return collect(Brng), ωrng, LDOS
end

function ldos_panel!(ax, name, side)
    Brng, ωrng, LDOS = load_ldos(name, side)
    heatmap!(ax, Brng, ωrng, LDOS; colormap = CMAP_LDOS, colorrange = LDOS_RANGE)
    hlines!(ax, [0]; color = LIGHT, linewidth = LW_SPINE)
    flux_lines!(ax; color = LIGHT)
    xlims!(ax, BLIMS...)
    return ax
end

"Secondary Φ/Φ₀ axis across the top of a panel."
function flux_axis!(gp)
    gax = Axis(gp; xaxisposition = :top,
               xlabel = rich(gk("Φ"), "/", gk("Φ"), sub("0")),
               backgroundcolor = :transparent)
    xlims!(gax, BLIMS...)
    gax.xticks = ([B_at_flux(n) for n in (0.5, 1, 1.5, 2)], ["0.5", "1", "1.5", "2"])
    hidespines!(gax)
    hideydecorations!(gax)
    hidexdecorations!(gax; ticks = false, ticklabels = false, label = false)
    return gax
end

function hero_fig3(; name = "C3")
    Trng, Ic, Is, Brng, Ict, _, Trng_Rs, Rs, Brng_Rs, Rs_theo, _, _, dVdImat =
        process_data(Exp_Data(), "fitted_$(name)")

    fig = Figure(size = HERO)
    gl_L = fig[1, 1] = GridLayout()
    gl_R = fig[1, 2] = GridLayout()

    # ---------------------------------------------------- (a) LDOS, left lead
    axA = Axis(gl_L[1, 1])
    ldos_panel!(axA, name, "left")
    ylims!(axA, 0, 0.3); axA.yticks = [0.2]
    hidexdecorations!(axA)
    text!(axA, 0.045, 0.035; text = "Left lead", color = LIGHT, align = (:left, :center),
          fontsize = FS_ANNOT, font = MONT_MED)
    for (edge, lbl) in zip(lobe_ends(), ("0L", "1L", "2L"))  # once per column,
        text!(axA, edge - 0.005, 0.265; text = lbl, color = LIGHT,   # in the empty
              align = (:right, :center), fontsize = FS_ANNOT, font = MONT_BOLD)  # top strip
    end
    flux_axis!(gl_L[1, 1])

    # --------------------------------------------------- (b) LDOS, right lead
    axB = Axis(gl_L[2, 1]; xlabel = "B (T)")
    ldos_panel!(axB, name, "right")
    ylims!(axB, -0.3, 0); axB.yticks = [-0.2]
    text!(axB, 0.045, -0.035; text = "Right lead", color = LIGHT, align = (:left, :center),
          fontsize = FS_ANNOT, font = MONT_MED)

    Label(gl_L[1:2, 0], rich(gk("ω"), " (meV)"); rotation = π/2,
          fontsize = FS_LABEL, color = GREY)
    Colorbar(gl_L[1:2, 2]; colormap = CMAP_LDOS, limits = (0, 1),
             ticks = ([0, 1], ["0", "1"]), label = "LDOS (a.u.)")

    # ------------------------------------------------------- (c) dV/dI and Ic
    axC = Axis(gl_R[1, 1]; ylabel = "I (nA)")
    hm = heatmap!(axC, Trng, Is, dVdImat; colormap = CMAP_C, colorrange = DVDI_RANGE,
                  colorscale = Makie.pseudolog10,
                  lowclip = cgrad(CMAP_C)[1], highclip = cgrad(CMAP_C)[end])
    # seismic is mostly light here, so the flux guides go dark rather than light
    flux_lines!(axC; color = GREY)
    # experiment first, theory over it — as the paper draws it, so the model curve
    # stays traceable through the marker band
    for s in (+1, -1)
        scatter!(axC, Trng, s .* Ic; color = C_EXP, markersize = MS_CURVE,
                 strokecolor = C_HALO, strokewidth = 1.5,
                 label = s == 1 ? "Experiment" : nothing)
    end
    for s in (+1, -1)
        lines!(axC, Brng, s .* Ict; color = C_HALO, linewidth = LW_HALO)
        lines!(axC, Brng, s .* Ict; color = C_THEO, linewidth = LW_CURVE,
               label = s == 1 ? "Theory" : nothing)
    end
    xlims!(axC, BLIMS...); ylims!(axC, -125, 125)
    axC.yticks = -100:50:100
    hidexdecorations!(axC; ticks = false)
    flux_axis!(gl_R[1, 1])

    # the argument: the Ic bump sits off the dashed lobe centre, and flips sign
    arrows2d!(axC, [0.118], [64], [0.212], [64]; argmode = :endpoint,
              color = C_ANNOT_ALT, shaftwidth = 4, tipwidth = 20, tiplength = 18)
    text!(axC, 0.165, 80; text = "skew", color = C_ANNOT_ALT, align = (:center, :bottom),
          fontsize = FS_ANNOT, font = MONT_BOLD)
    arrows2d!(axC, [0.370], [44], [0.276], [44]; argmode = :endpoint,
              color = C_ANNOT, shaftwidth = 4, tipwidth = 20, tiplength = 18)
    text!(axC, 0.323, 60; text = "skew", color = C_ANNOT, align = (:center, :bottom),
          fontsize = FS_ANNOT, font = MONT_BOLD)

    # Makie draws the legend backdrop as part of the frame, so the frame stays
    # on: grey text on a mid-tone heatmap is otherwise unreadable at 2 m.
    axislegend(axC; position = :rb, orientation = :horizontal, colgap = 22,
               backgroundcolor = LIGHT, framevisible = true,
               framecolor = RULE, framewidth = LW_SPINE, padding = (14, 14, 6, 6))
    Colorbar(gl_R[1, 2], hm; label = rich("dV/dI (k", gk("Ω"), ")"),
             ticks = (collect(DVDI_RANGE), ["0.5", "1.5"]))

    # ---------------------------------------------------------- (d) R₀ / R_N
    axD = Axis(gl_R[2, 1]; xlabel = "B (T)",
               ylabel = rich("R", sub("0"), "/R", sub("N")))
    flux_lines!(axD; color = GREY)
    scatter!(axD, Trng_Rs, Rs; color = C_EXP, markersize = MS_CURVE,
             strokecolor = C_HALO, strokewidth = 1.5)
    lines!(axD, Brng_Rs, Rs_theo; color = C_HALO, linewidth = LW_HALO)
    lines!(axD, Brng_Rs, Rs_theo; color = C_THEO, linewidth = LW_CURVE)
    xlims!(axD, BLIMS...); ylims!(axD, -0.04, 1.14)
    axD.yticks = [0, 0.5, 1]
    for (edge, lbl) in zip(lobe_ends(), ("0L", "1L", "2L"))
        text!(axD, edge - 0.005, 0.1; text = lbl, color = GREY, align = (:right, :center),
              fontsize = FS_ANNOT, font = MONT_BOLD)
    end

    # ------------------------------------------------------------- geometry --
    rowgap!(gl_L, 1, 10)                   # a clear white split between (a) and (b)
    rowsize!(gl_R, 1, Relative(0.60))
    rowgap!(gl_R, 1, 14)
    colgap!(gl_L, 1, 4); colgap!(gl_L, 2, 12)
    colgap!(gl_R, 1, 12)
    colgap!(fig.layout, 1, 46)
    colsize!(fig.layout, 1, Relative(0.5))

    panel_label!(gl_L[1, 1, TopLeft()], "a")
    panel_label!(gl_L[2, 1, TopLeft()], "b")
    panel_label!(gl_R[1, 1, TopLeft()], "c")
    panel_label!(gl_R[2, 1, TopLeft()], "d")

    return fig
end

"""
    exact_page!(path, size_pt)

CairoMakie rounds the scene to whole points, so a 771 × 240 mm slot comes out as
a 2186 × 680 pt page — 0.17 mm too wide, 0.11 mm too short. Reset the MediaBox to
the exact slot size. qpdf's QDF form is edited in place and repaired with fix-qdf,
so every content stream survives byte for byte; the only thing the width loses is
0.5 pt of empty figure padding.
"""
function exact_page!(path::AbstractString, size_pt::Tuple)
    qpdf, fixqdf = Sys.which("qpdf"), Sys.which("fix-qdf")
    if qpdf === nothing || fixqdf === nothing
        @warn "qpdf/fix-qdf not on PATH — page stays at CairoMakie's rounded size"
        return false
    end
    mktempdir() do dir
        qdf = joinpath(dir, "qdf.pdf")
        run(pipeline(`$qpdf --qdf --object-streams=disable $path $qdf`; stdout = devnull))
        bytes = read(qdf)
        i = first(something(findfirst(codeunits("/MediaBox"), bytes)))
        j = first(something(findnext(codeunits("]"), bytes, i)))
        box = "/MediaBox [ 0 0 $(round(size_pt[1], digits = 4)) $(round(size_pt[2], digits = 4)) ]"
        write(qdf, vcat(bytes[1:i-1], codeunits(box), bytes[j+1:end]))
        fixed = joinpath(dir, "fixed.pdf")
        open(fixed, "w") do io
            run(pipeline(`$fixqdf $qdf`; stdout = io))
        end
        # QDF form is uncompressed — put the streams back the way Cairo had them
        run(pipeline(`$qpdf --stream-data=compress --object-streams=generate $fixed $path`;
                     stdout = devnull))
    end
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    fig = hero_fig3()
    mkpath("poster/figures")
    save("poster/figures/hero_fig3.pdf", fig; pt_per_unit = 1)
    save("poster/figures/hero_fig3.png", fig; px_per_unit = 150/72)   # preview only
    exact_page!("poster/figures/hero_fig3.pdf", HERO)
    @info "wrote poster/figures/hero_fig3.{pdf,png}"
end
