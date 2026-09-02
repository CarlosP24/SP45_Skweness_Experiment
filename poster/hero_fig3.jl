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

# --------------------------------------------------------------------------
# §5 additions, drawn only when `labelled = true`. Everything else about the
# hero is untouched, and the unlabelled export is unchanged.
#
# The hC/eC identification is NOT eyeballed off the map. It comes from the
# per-m_J LDOS of the same wire (single-mode `build_cyl`, one sector per run),
# which puts, at Φ = 0.90 in the right lead:
#
#     m_J = ±1  →  ω = −0.136 meV   hole-like     (hC throughout the 1L lobe)
#     m_J = ±3  →  ω = −0.040 meV   electron-like (eC until it crosses zero)
#
# The ±3 analog is the one that crosses zero energy at Φ ≈ 1.05–1.075 and
# converts to hC — the same conversion the CPR figure's panels (b) and (c)
# bracket. So the eC arrow here points at exactly the analog that figure is about.
#
# Caveat worth knowing: the map itself is the multimode model summed over m_J,
# while the branch identification is single-mode per sector. The two were
# validated against each other (crossing fluxes agree for all four ±m_J pairs);
# see data/theory/PROVENANCE.md.
const Φ0_PER_T = 6.2229          # Φ/Φ₀ per tesla, for R_LP = 64 nm
B_of_Φ(Φ) = Φ / Φ0_PER_T

const LBL_hC_Φ    = 0.95            # flux at which the two analogs are labelled
const LBL_hC_ω = 0.14
const LBL_eC_Φ    = 1.05            # flux at which the two analogs are labelled
const LBL_eC_ω = 0.09

# The three fluxes the CPR mechanism figure shows, with its symbols (§5).
const MARK_Φ   = (0.60, 1.05, 1.10)
const MARK_SYM = (:circle, :rect, :utriangle)

"""
    flux_marks!(ax, y)

The three fluxes the CPR mechanism figure shows, at height `y` in `ax`'s own
units, using that figure's symbols.

Light fill with a dark stroke, not the reverse: the row crosses both the pale
band near the 0L/1L edge in (c) and the dark ground of the LDOS maps, and this
way round reads on both. The square and triangle sit 0.008 T apart because those
two panels deliberately bracket one zero-energy crossing.
"""
function flux_marks!(ax, y)
    for (Φ, sym) in zip(MARK_Φ, MARK_SYM)
        scatter!(ax, [B_of_Φ(Φ)], [y]; marker = sym, markersize = 20,
                 color = LIGHT, strokecolor = C_ANNOT, strokewidth = 2.5)
    end
    return ax
end

"""
Short leader arrow plus label, in the hC/eC colour code.

`side` says which way the label sits from the feature. The eC label has to come
in from the right: at ω ≈ −0.04 a left-hand label runs straight into the
existing "Right lead" annotation, which sits at the same height.
"""
function analog_label!(ax, Bfeat, ωfeat, text_str, colour; side = :left, dB = 0.032)
    s  = side === :left ? -1 : +1
    B0 = Bfeat + s * dB
    arrows2d!(ax, [B0 - s * 0.008], [ωfeat], [Bfeat + s * 0.004], [ωfeat];
              argmode = :endpoint, color = colour,
              shaftwidth = 3.5, tipwidth = 16, tiplength = 14)
    text!(ax, B0 - s*0.004, ωfeat; text = text_str, color = colour,
          align = (side === :left ? :right : :left, :center),
          fontsize = FS_ANNOT, font = MONT_BOLD)
    return ax
end

function hero_fig3(; name = "C3", labelled = false)
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

    if labelled
        # Both labels go in (b) rather than one per panel: the two families are
        # only distinguishable in the negative-energy sector, and (a) shows ω > 0.
        analog_label!(axA, B_of_Φ(LBL_eC_Φ), LBL_eC_ω, "eC", C_EC; side = :right)
        analog_label!(axA, B_of_Φ(LBL_hC_Φ), LBL_hC_ω, "hC", C_HC; side = :left)
        # Repeated here as well as in (c), so the three CPR panels are tied to the
        # spectrum they come from and not only to the measured critical current.
        # Below the gap edge, where the map is dark and empty.
        flux_marks!(axB, -0.285)
    end

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

    # Flux markers tying (c) to the CPR figure's three panels, matching (b)'s row
    labelled && flux_marks!(axC, -108)

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

if abspath(PROGRAM_FILE) == @__FILE__
    mkpath("poster/figures")
    # Named for the paper figure it rebuilds, per this repo's convention. The
    # build doc calls this hero_fig1; poster/CLAUDE.md, MANIFEST.md and the file
    # already on disk all call the hero slot fig3, so fig3 is kept.
    for (suffix, labelled) in (("", false), ("_labelled", true))
        fig = hero_fig3(; labelled)
        pdf = "poster/figures/hero_fig3$(suffix).pdf"
        save(pdf, fig; pt_per_unit = 1)
        save("poster/figures/hero_fig3$(suffix).png", fig; px_per_unit = 150/72)  # preview only
        exact_page!(pdf, HERO)
        @info "wrote poster/figures/hero_fig3$(suffix).{pdf,png}"
    end
end
