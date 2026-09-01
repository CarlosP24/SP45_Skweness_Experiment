# Second supporting slot — 298.4 × 112.5 mm (80 % of the supporting slot's width,
# 75 % of its height). Paper Fig. 2(b): dV/dI against axial field and temperature,
# where the core's critical temperature carries the same skew as the critical current.
#
# Re-style only. The map is the paper's own `prepare_data` output, unmodified.
#
# Run from anywhere:  julia --project=plots poster/support2_fig2b.jl

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)   # the paper's readers use paths relative to the repo root

using Parameters, DelimitedFiles, SavitzkyGolay, Statistics
using CairoMakie

include(joinpath(ROOT, "poster", "style", "poster_theme.jl"))
include(joinpath(ROOT, "plots", "kitchen", "data_reader.jl"))   # paper code, unmodified

set_theme!(poster_theme)

const BLIMS      = (0.0, 0.4)        # T — as in the paper, not one gauss wider
const DVDI_RANGE = (4e2, 1.5e3)      # Ω
const CMAP_C     = :seismic          # the paper's own dV/dI map, kept as published

# T = 250 mK sweep, V_G = 80 V — the same device data as paper Fig. 2(b)
const EXP_T = Exp_Data(; keysight = "data_exp/250mK_DEV13_294_35_V_drop.dat",
                         lockin   = "data_exp/250mK_DEV13_294_35_R_demod.dat")

function support2_fig2b()
    # first column of this file is T, second is B — hence the (Ts, Bs) order
    Ts, Bs, dVdImat, _, _ = prepare_data(EXP_T, 0.6; Icoef = 1)
    iB = findmin(abs.(Bs .- BLIMS[2]))[2]
    dVdI = abs.(dVdImat[:, 1:iB])'

    fig = Figure(size = SUPPORT2)
    ax = Axis(fig[1, 1]; xlabel = "B (T)", ylabel = "T (K)")
    hm = heatmap!(ax, Bs[1:iB], Ts, dVdI; colormap = CMAP_C, colorrange = DVDI_RANGE,
                  lowclip = cgrad(CMAP_C)[1], highclip = cgrad(CMAP_C)[end])
    # this map is saturated dark at both ends, so the flux guides go light
    flux_lines!(ax; color = LIGHT)
    xlims!(ax, BLIMS...)
    ax.yticks = [0.5, 1.0]

    for (edge, lbl) in zip(lobe_ends(), ("0L", "1L", "2L"))   # in the normal-state strip
        text!(ax, edge - 0.005, 1.36; text = lbl, color = LIGHT,
              align = (:right, :center), fontsize = FS_ANNOT, font = MONT_BOLD)
    end

    # the same claim the hero makes, in temperature: the surviving core lobe sits
    # right of 1L's centre, not on it
    arrows2d!(ax, [0.118], [0.72], [0.207], [0.72]; argmode = :endpoint,
              color = C_ANNOT, shaftwidth = 4, tipwidth = 20, tiplength = 18)
    text!(ax, 0.162, 0.80; text = "skew", color = C_ANNOT, align = (:center, :bottom),
          fontsize = FS_ANNOT, font = MONT_BOLD)

    Colorbar(fig[1, 2], hm; label = rich("dV/dI (k", gk("Ω"), ")"),
             ticks = (collect(DVDI_RANGE), ["0.4", "1.5"]))
    colgap!(fig.layout, 1, 12)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    fig = support2_fig2b()
    mkpath("poster/figures")
    save("poster/figures/support2_fig2b.pdf", fig; pt_per_unit = 1)
    save("poster/figures/support2_fig2b.png", fig; px_per_unit = 150/72)   # preview only
    exact_page!("poster/figures/support2_fig2b.pdf", SUPPORT2)
    @info "wrote poster/figures/support2_fig2b.{pdf,png}"
end
