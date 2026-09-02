# CPR mechanism figure — 373 × 143 mm column slot (build doc §3 slot).
#
# What this figure says, and why it exists: the hero shows that Ic is skewed and
# that theory reproduces it. It does not show WHY a skewed spectrum skews the
# supercurrent. This does.
#
# The mechanism, as the calculation actually has it (see figures/diagnostics/):
# within the 1L lobe, electron-like CdGM analogs (eC, π-junction-like) convert to
# hole-like ones (hC, 0-junction-like) in symmetric ±m_J pairs, innermost first,
# at four discrete fluxes. Each conversion removes a negative contribution and
# adds a positive one, so the total steps UP; those upward steps push the Ic
# maximum to Φ = 1.0875, past the lobe centre at 1.0 — the rightward skew.
#
# NOTE: this is the opposite direction to the build doc's §1 and §9 draft, which
# describe hC → eC and the total stepping down. The data is the authority here;
# the direction was checked numerically (ΔIc = +0.079, +0.156, +0.077 at the
# three crossings inside the lobe).
#
# Run:  julia --project=plots poster/cpr_mechanism.jl

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using Quantica, FullShell, JLD2, Parameters, CairoMakie
using Printf

include(joinpath(ROOT, "src", "utilities.jl"))          # filter_πs
include(joinpath(ROOT, "src", "models", "params.jl"))
include(joinpath(ROOT, "src", "models", "wires.jl"))
include(joinpath(ROOT, "src", "models", "junctions.jl"))
include(joinpath(ROOT, "src", "models", "cpr.jl"))
include(joinpath(ROOT, "poster", "style", "poster_theme.jl"))

set_theme!(poster_theme)

const CPR_PATH = "data/theory/cpr/cpr_fitted_C3.jld2"

# The three fluxes, fixed by the §3 diagnostic and not before it: two bracketing
# the ±3 conversion at Φ = 1.075 (the largest step, ΔIc = +0.156) and one well
# away from any crossing. Ordered left to right in increasing flux, matching the
# hero's field axis direction.
const PANEL_Φ   = (0.60, 1.05, 1.10)
const PANEL_SYM = (:circle, :rect, :utriangle)   # repeated on the hero, §5
const PANEL_LET = ("a", "b", "c")

# Well above the 0.5 mm floor (0.5 mm = 1.42 pt). With the total gone there is no
# longer a thick curve for these to be read against, so they carry the figure and
# can afford the weight.
const LW_SECTOR = 5.0     # 1.76 mm

# Supercurrent in nA, on the hero's own conversion so the two figures can be read
# against each other: `ics()` in plots/plotters/theo_exp_core.jl scales the raw
# Josephson current by conv = 74.4 nA, and process_data then applies pref = 0.8.
#
# The hero's theory curve is at 460 mK and this is at T = 0, so these values run
# ~20 % above the hero's at the same flux (the diagnostic put the ratio at 0.837).
# That is the deliberate temperature difference, not a units mismatch.
const NA_PER_UNIT = 74.4 * 0.8

"Index of the flux in `Φrng` closest to `Φ`."
nearest(Φrng, Φ) = argmin(abs.(Φrng .- Φ))

"""
SP43's classification (plots/plotters/plot_cphase.jl:33-35): a sector is
hole-like if its current peaks at φ < π (0-junction-like) and electron-like if
it peaks at φ > π (π-junction-like).

Reused rather than reimplemented, per build doc §2. Cross-checked against the
sign-of-Ic criterion §2 describes: the two agree on 99.62 % of sector-flux
points, disagreeing only where a sector sits exactly on a boundary.
"""
function classify(Jz, φrng)
    all(isnan, Jz) && return :nan
    return φrng[argmax(replace(Jz, NaN => -Inf))] < π ? :hC : :eC
end

"""
    cpr_mechanism()

Three panels, one shared symmetric scale, individual m_J contributions only.

The total is deliberately absent. It reaches 26 nA while the sector that actually
distinguishes panel (b) from panel (c) — m_J = ±3, converting from electron-like
to hole-like — is 2.1 nA on one side of the crossing and 2.9 nA on the other. Any
linear axis wide enough for the total flattens that flip to a tenth of the
half-height, where it cannot be read at three metres. Dropping the total lets the
largest CdGM contribution set the scale instead, and the hC/eC counts become
countable, which is this figure's whole job. The critical current itself is
already on the hero, as the theory curve in panel (c).
"""
function cpr_mechanism()
    res  = load(CPR_PATH)["res"]
    sys  = res.params
    Js   = res.Js
    Zs   = sort(collect(keys(Js)))
    Φrng = collect(sys.Φrng)
    φrng = filter_πs(collect(sys.φrng))

    # Scale set by the largest individual CdGM contribution across all three
    # panels — m_J = 0 at Φ = 1.05, 8.9 nA — not by any total.
    gmax = 0.0
    for Φt in PANEL_Φ
        i = nearest(Φrng, Φt)
        for Z in Zs
            v = filter(isfinite, Js[Z][i, :])
            isempty(v) || (gmax = max(gmax, maximum(abs, v)))
        end
    end
    gmax *= NA_PER_UNIT
    ylim = 1.08 * gmax

    fig = Figure(size = CPR)
    gl  = fig[1, 1] = GridLayout()

    for (k, Φt) in enumerate(PANEL_Φ)
        i  = nearest(Φrng, Φt)
        Φ  = Φrng[i]
        # π has no Montserrat glyph in any release, so every one of them — tick
        # labels included — goes through gk(). Left as a plain string it renders,
        # but silently drags TeX Gyre Heros into the PDF as a fifth font.
        ax = Axis(gl[1, k];
                  xlabel = k == 2 ? rich(gk("φ")) : "",
                  ylabel = k == 1 ? rich("J", sub("S"), " (nA)") : "",
                  xticks = ([0, π, 2π],
                            [rich("0"), rich(gk("π")), rich("2", gk("π"))]))

        hlines!(ax, [0]; color = (GREY, 0.55), linewidth = LW_SPINE)

        # individual sectors, thin; hC and eC drawn in the poster-wide code
        nhC = neC = 0
        for Z in Zs
            Jz = Js[Z][i, :] .* NA_PER_UNIT
            c  = classify(Jz, φrng)
            c === :nan && continue
            c === :hC ? (nhC += 1) : (neC += 1)
            lines!(ax, φrng, Jz;
                   color = c === :hC ? C_HC : C_EC, linewidth = LW_SECTOR)
        end

        # The critical current of the *total*, which this figure no longer draws —
        # stated so the step across the ±3 conversion stays on the page even though
        # the curve carrying it does not.
        Ic = maximum(replace(sum(Js[Z][i, :] for Z in Zs) .* NA_PER_UNIT,
                             NaN => -Inf))

        xlims!(ax, 0, 2π)
        ylims!(ax, -ylim, ylim)

        # Flux value above the panel, not inside it — the panel interior is
        # entirely given over to the curves now.
        Label(gl[1, k, Top()],
              rich(gk("Φ"), "/", gk("Φ"), sub("0"), " = $(round(Φ; digits = 2))");
              fontsize = FS_ANNOT, font = MONT_MED, color = GREY,
              padding = (0, 0, 6, 0))

        # the symbol that ties this panel to its flux markers on the hero —
        # drawn identically: LIGHT fill, C_ANNOT stroke, same size and weight
        scatter!(ax, [0.945], [0.93]; space = :relative, marker = PANEL_SYM[k],
                 markersize = 20, color = LIGHT,
                 strokecolor = C_ANNOT, strokewidth = 2.5)

        # Counts and Ic go in the lower-left. Between φ = 0 and π every hC curve is
        # positive and the eC ones only reach about −2 nA, so the band below −6 nA
        # is empty across more than half the panel — the only corner with room now
        # that the sectors fill the frame. Counts carry their family colour, so the
        # code is legible without reference to the legend.
        text!(ax, 0.05, 0.175;
              text = rich(rich("$nhC hC"; color = C_HC), "   ",
                          rich("$neC eC"; color = C_EC)),
              space = :relative, align = (:left, :bottom),
              fontsize = FS_LEGEND, font = MONT_MED)
        text!(ax, 0.05, 0.045;
              text = rich("I", sub("c"), " = ", @sprintf("%.1f", Ic), " nA"),
              space = :relative, align = (:left, :bottom),
              fontsize = FS_LEGEND, font = MONT_REG, color = GREY)

        # One scale across all three panels, so the numbers mean something and the
        # left panel carries the ticks (§4).
        #
        # Explicit labels with an ASCII hyphen: Makie's automatic numeric ticks use
        # a Unicode minus (U+2212), which Montserrat has no glyph for, so leaving
        # them automatic silently drags TeX Gyre Heros into the PDF as a fifth
        # font — the same trap as the π tick labels above.
        showticks = k == 1
        showticks && (ax.yticks = ([-5, 0, 5], [rich("-5"), rich("0"), rich("5")]))
        hideydecorations!(ax; label = k != 1,
                          ticks = !showticks, ticklabels = !showticks,
                          grid = true, minorgrid = true, minorticks = true)
        panel_label!(gl[1, k, TopLeft()], PANEL_LET[k]; padding = (0, 10, -30, 0))
    end

    # hC / eC label once, not repeated across panels (§4)
    Legend(gl[1, 1, Bottom()],
           [LineElement(color = C_HC, linewidth = LW_SECTOR),
            LineElement(color = C_EC, linewidth = LW_SECTOR)],
           ["hC", "eC"];
           orientation = :horizontal, framevisible = false,
           labelsize = FS_LEGEND, labelcolor = GREY,
           padding = (0, 0, 0, 0), colgap = 20, patchsize = (40, 18),
           tellheight = false, tellwidth = false,
           halign = :center, valign = :bottom,
           margin = (0, 0, 0, -34))

    colgap!(gl, mm2pt(4))
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    # The build doc §8 says figures/poster/. This repo already keeps every poster
    # deliverable in poster/figures/, next to MANIFEST.md, and having both
    # figures/poster/ and poster/figures/ is a trap nobody needs. Repo convention
    # wins; the §3 diagnostic still goes to figures/diagnostics/ as the doc asks,
    # since that one is not a poster deliverable.
    mkpath("poster/figures")
    fig = cpr_mechanism()
    save("poster/figures/cpr_mechanism.pdf", fig; pt_per_unit = 1)
    save("poster/figures/cpr_mechanism.png", fig; px_per_unit = 300/72)
    exact_page!("poster/figures/cpr_mechanism.pdf", CPR)
    @info "wrote poster/figures/cpr_mechanism.{pdf,png}"
end
