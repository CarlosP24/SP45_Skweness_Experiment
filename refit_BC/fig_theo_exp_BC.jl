# Step 6 of the devices-B/C fit: Fig. 3 of the preprint (plots/fig_theo_exp.jl) for devices B and C.
#   julia --project=plots refit_BC/fig_theo_exp_BC.jl [T_e mode] [rank] [grid names...]   (default grids B C)
# Same layout and styling as fig_theo_exp: (a, b) LDOS of the left and right leads (positive / negative
# energies), (c) low-current dV/dI with I_c^AH (green) and I_c^th (orange), (d) R_0/R_N, data (black) and
# the Ambegaokar–Halperin curve of I_c^th (purple). Differences: the theory is the m_J-resolved sum with
# the physical Zeeman (refit_BC/final.jl) and the data are the digitized SI maps (refit_BC/target.jl).
# Output: refit_BC/figures/fig_theo_exp_<grid>.pdf (+ .png)
using Parameters, JLD2, CairoMakie, Statistics, Bessels, Printf, DelimitedFiles, SavitzkyGolay, Roots, Peaks
include(joinpath(pwd(), "plots", "kitchen", "data_reader.jl"))
include(joinpath(pwd(), "plots", "kitchen", "Ic_reader.jl"))            # γAH
include(joinpath(pwd(), "plots", "plotters", "theo_exp_core.jl"))      # panel_dVdI, panel_resistance
include(joinpath(pwd(), "refit_BC", "devices.jl"))

const MODE = isempty(ARGS) ? "fixed" : ARGS[1]
const RANK = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
const GRIDS = length(ARGS) >= 3 ? ARGS[3:end] : ["B", "C"]
const CONV = 74.4

function panel_LDOS_sm(ax, Bs, ωs, LDOS, side; colorrange)
    ω = real.(ωs)
    L = copy(LDOS)
    if side == "left"                       # positive energies for the left lead, as panel_LDOS
        L = reverse(L; dims = 2); ω = -reverse(ω)
    end
    heatmap!(ax, Bs, ω, L; colormap = :thermal, colorrange, rasterize = 5)
    hlines!(ax, [0]; color = :white, linewidth = 2)
end

function fig_theo_exp_dev(name; mode = MODE, rank = RANK)
    fn = load(joinpath("refit_BC", "out", "final_$(name)_$(mode)_$(rank).jld2"))
    dev = fn["dev"]
    tg = load(joinpath("refit_BC", "out", "target_$dev.jld2"))
    G = GEOMETRY[dev]
    R = G.R + G.d / 2                        # R_LP, for the flux axis (Φ = π R_LP² B / Φ0, as Fig. 3)
    Φ0 = 2067.833
    Bmax = maximum(fn["Bs"])
    nmax = floor(Int, Bmax * π * R^2 / Φ0 + 0.5)
    Bns = [n * Φ0 / (π * R^2) for n in 1:nmax if n * Φ0 / (π * R^2) < Bmax]
    Bhalfs = [(n + 0.5) * Φ0 / (π * R^2) for n in 0:nmax if (n + 0.5) * Φ0 / (π * R^2) < Bmax]
    # scale and T_e from the current scoring of the same grid (analyze.jl), matched by candidate
    sc = load(joinpath("refit_BC", "out", "$(name)_$(mode)_scores.jld2"))
    ic = findfirst(==(fn["cand"]), sc["cands"])
    Te, s = ic === nothing ? (fn["Te"], fn["s"]) : (sc["Te"][ic], sc["s"][ic])

    # theory: I_c^th in nA with the fitted scale, and its AH resistance
    Bth = fn["Bs"]
    Ict = CONV * s .* fn["Ic"]
    Rs_theo = @. besseli0(Ict / (0.0209 * Te) / 2)^(-2)

    # experiment: lobes n ≥ 1 only, as Fig. 3 (no R_N in 0L)
    B, r0 = tg["B"], tg["r0"]
    iB = findall(B .<= Bmax)
    m1 = [b > tg["edges"][1] for b in B[iB]]
    IcAH = 0.0209 * Te .* γAH.(abs.(r0[iB]))
    Ic_plot = [m ? x : NaN for (m, x) in zip(m1, IcAH)]

    fig = Figure(size = (600, 750), fontsize = 20)
    # LDOS colour range: same relative scale for both leads (published Fig. 3 uses (0, 2e-3))
    cr = (0, quantile(filter(isfinite, vec(vcat(fn["LDOS_L"], fn["LDOS_R"]))), 0.97))
    for (row, side, key) in ((1, "left", "LDOS_L"), (2, "right", "LDOS_R"))
        ax = Axis(fig[row, 1];)
        panel_LDOS_sm(ax, fn["BLDOS"], fn["ωs"], fn[key], side; colorrange = cr)
        hidexdecorations!(ax, ticks = false, grid = false)
        hideydecorations!(ax)
        vlines!(ax, Bns; color = :white, linestyle = :dash)
        vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)
        xlims!(ax, 0, Bmax)
        if row == 1
            for n in 0:nmax
                x = (n + 0.4) * Φ0 / (π * R^2)
                # skip 0L where it would run into the "Left" label (narrow lobes of device C)
                (n == 0 && x < 0.15 * Bmax) && continue
                x < Bmax && text!(ax, x, 0.25; text = "$(n)L", align = (:center, :center), color = :white)
            end
            text!(ax, 0.01 * Bmax / 0.4, 0.25; text = "Left", color = :white, align = (:left, :center), fontsize = 20)
            gax = Axis(fig[1, 1]; xaxisposition = :top, xlabel = L"\Phi / \Phi_0")
            xlims!(gax, 0, Bmax)
            ticks = sort(vcat([(b, n) for (n, b) in enumerate(Bns)], [(b, n - 0.5) for (n, b) in enumerate(Bhalfs)]); by = first)
            gax.xticks = (first.(ticks), [x == round(x) ? L"%$(Int(x))" : L"%$(x)" for x in last.(ticks)])
            hidespines!(gax); hideydecorations!(gax)
            hidexdecorations!(gax, ticks = false, ticklabels = false, label = false)
        else
            text!(ax, 0.01 * Bmax / 0.4, -0.25; text = "Right", color = :white, align = (:left, :center), fontsize = 20)
        end
    end
    ax = Axis(fig[1:2, 1]; ylabel = L"$\omega$ (meV)", ylabelpadding = -35)
    hidexdecorations!(ax); ylims!(ax, (-0.3, 0.3)); ax.yticks = [-0.2, 0.2]; hidespines!(ax)
    Colorbar(fig[1:2, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    colorrange = (5e2, 1.5e3)
    ax = Axis(fig[3, 1]; xlabel = L"$B$ (T)", ylabel = L"$I$ (nA)")
    lims, hmap = panel_dVdI(ax, tg["dVdImat"][iB, :], B[iB], tg["Is"], Ic_plot, Bth, Ict; colormap = :seismic,
        lowclip = cgrad(:seismic)[1], colorrange, highclip = cgrad(:seismic)[end])
    ax.ylabelpadding = -18
    xlims!(ax, 0, Bmax)
    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)
    Imax = maximum(tg["Is"])
    for n in 0:nmax
        x = (n + 0.3) * Φ0 / (π * R^2)
        x < Bmax && text!(ax, x, -0.9 * Imax; text = "$(n)L", align = (:center, :center), color = :white)
    end
    hidexdecorations!(ax, ticks = false, grid = false)
    Colorbar(fig[3, 2], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20,
        ticks = ([first(colorrange), last(colorrange)], [first(colorrange), last(colorrange)] ./ 1e3 .|> string))

    ax = Axis(fig[4, 1]; xlabel = L"$B$ (T)", ylabel = L"$R_0/R_N$")
    panel_resistance(ax, B[iB][m1], r0[iB][m1], Bth, Rs_theo; legend = false)
    # legend in the empty corner of each device (Fig. 3 uses (1, 0.5), where these devices have data)
    if dev == "C"     # compact, in the empty top-left corner (0L), above the 1L-edge spike of the theory
        axislegend(ax; position = (0.0, 1.0), labelsize = 14, patchsize = (14, 10), padding = (5, 5, 2, 2), rowgap = 0)
    else
        axislegend(ax; position = (0.62, 0.06))
    end
    xlims!(ax, 0, Bmax)
    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)
    for n in 0:nmax
        x = (n + 0.4) * Φ0 / (π * R^2)
        x < Bmax && text!(ax, x, 0.1; text = "$(n)L", align = (:center, :center), color = :black)
    end

    rowsize!(fig.layout, 1, Relative(0.25 / 2)); rowsize!(fig.layout, 2, Relative(0.25 / 2)); rowsize!(fig.layout, 4, Relative(0.25))
    colgap!(fig.layout, 1, 5); rowgap!(fig.layout, 1, 5); rowgap!(fig.layout, 2, 5); rowgap!(fig.layout, 3, 5)
    style = (font = "CMU Serif Bold", fontsize = 20)
    for (r, l) in enumerate(("a", "b", "c", "d"))
        Label(fig[r, 1, TopLeft()], l, padding = (-40, 0, -30, 0); style...)
    end
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    for name in GRIDS
        isfile(joinpath("refit_BC", "out", "final_$(name)_$(MODE)_$(RANK).jld2")) || (println("no final curves for $name yet"); continue)
        fig = fig_theo_exp_dev(name)
        save(joinpath("refit_BC", "figures", "fig_theo_exp_$name.pdf"), fig)
        save(joinpath("refit_BC", "figures", "fig_theo_exp_$name.png"), fig; px_per_unit = 2)
        println("saved refit_BC/figures/fig_theo_exp_$name.{pdf,png}")
    end
end
