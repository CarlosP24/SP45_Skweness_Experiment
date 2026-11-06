# §3 diagnostic pass for the CPR mechanism figure.
#
# NOT a poster figure. This decides two things before any figure is designed:
#
#   1. Does the single-mode m_J-resolved model reproduce the hero? (Σ_Z Ic vs the
#      hero's cached Ic.) If not, a convention was missed and nothing downstream
#      is trustworthy.
#
#   2. Is the 1L lobe a COMPETITION story or a SUBSTITUTION one? That choice
#      determines which three fluxes the figure's panels show, and the build doc
#      forbids picking them before this runs.
#
# Run:  julia --project=plots plots/cpr_diagnostic.jl

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using Quantica, FullShell, JLD2, Parameters, CairoMakie

"""
Linear interpolation of `(xs, ys)` at `x`, clamped at the ends.

Hand-rolled rather than pulling `Interpolations` into `plots/Project.toml`: that
environment is what builds the hero, and it does not need a new dependency for
five lines in a diagnostic. `xs` must be sorted ascending.
"""
function lininterp(xs, ys, x)
    x <= first(xs) && return first(ys)
    x >= last(xs)  && return last(ys)
    j = searchsortedlast(xs, x)
    j >= length(xs) && return last(ys)
    t = (x - xs[j]) / (xs[j+1] - xs[j])
    return (1 - t) * ys[j] + t * ys[j+1]
end

# The .jld2 holds CPR_System and Params objects; JLD2 needs the definitions.
include(joinpath(ROOT, "src", "utilities.jl"))
include(joinpath(ROOT, "src", "models", "params.jl"))
include(joinpath(ROOT, "src", "models", "wires.jl"))
include(joinpath(ROOT, "src", "models", "junctions.jl"))
include(joinpath(ROOT, "src", "models", "cpr.jl"))

const CPR_PATH  = "data/theory/cpr/cpr_fitted_C3.jld2"
const HERO_PATH = "data/Js/fitted_C3.jld2"
const OUTDIR    = "figures/diagnostics"

# SP43's own colours (plots/fig_jos_triv.jl): turquoise = 0-junction-like (hC),
# coral = π-junction-like (eC). Kept verbatim here — this is a diagnostic, not
# the poster, so the coral/orange clash does not apply.
const TURQ  = RGBf(103/255, 203/255, 194/255)
const CORAL = RGBf(255/255, 127/255,  85/255)

res = load(CPR_PATH)["res"]
sys = res.params
Js  = res.Js
Zs  = sort(collect(keys(Js)))

Φrng = collect(sys.Φrng)
deleteat!(Φrng, findall(x -> isapprox(x, 0), Φrng))
φrng = filter_πs(collect(sys.φrng))

nΦ, nφ = size(Js[first(Zs)])
@info "loaded" nΦ nφ Zs

# --------------------------------------------------------------------------
# 1. Classification, using SP43's actual criterion
#
# SP43 (plots/plotters/plot_cphase.jl:33-35) classifies a sector by WHERE its
# current peaks in phase, not by the sign of its Ic:
#
#     φM = φrng[argmax(J_Z)] ;  icol = ceil(Int, φM/π)     # 1 → hC, 2 → eC
#
# For a near-sinusoidal CPR that agrees with the sign convention the build doc
# describes; for a skewed one it need not. Reusing SP43's, as §2 instructs.
# --------------------------------------------------------------------------
"Return :hC or :eC for sector `Z` at flux index `i`, by SP43's rule."
function classify(Z, i)
    Jz = Js[Z][i, :]
    all(isnan, Jz) && return :nan
    φM = φrng[argmax(replace(Jz, NaN => -Inf))]
    return φM < π ? :hC : :eC
end

# --------------------------------------------------------------------------
# 2. Weights
#
# Two different things are worth plotting and they answer different questions.
#
#   * peak weight   max_φ J_Z(φ) — how much current the sector could carry.
#     Always positive, so it cannot show cancellation.
#
#   * signed weight J_Z(φ*), where φ* maximises the TOTAL. This is the sector's
#     actual contribution to Ic, and it is negative for π-junction sectors. This
#     is the one that shows an hC being lost and an eC taking its place.
# --------------------------------------------------------------------------
Jtot = sum(Js[Z] for Z in Zs)
Ic   = [all(isnan, Jtot[i, :]) ? NaN : maximum(replace(Jtot[i, :], NaN => -Inf)) for i in 1:nΦ]
φstar = [all(isnan, Jtot[i, :]) ? 1 : argmax(replace(Jtot[i, :], NaN => -Inf)) for i in 1:nΦ]

hC_peak = zeros(nΦ); eC_peak = zeros(nΦ)
hC_sgn  = zeros(nΦ); eC_sgn  = zeros(nΦ)
for i in 1:nΦ, Z in Zs
    c = classify(Z, i)
    c === :nan && continue
    Jz = Js[Z][i, :]
    pk = maximum(replace(Jz, NaN => -Inf))
    sg = Jz[φstar[i]]
    isfinite(pk) || continue
    if c === :hC
        hC_peak[i] += pk
        isfinite(sg) && (hC_sgn[i] += sg)
    else
        eC_peak[i] += pk
        isfinite(sg) && (eC_sgn[i] += sg)
    end
end

nhC = [count(Z -> classify(Z, i) === :hC, Zs) for i in 1:nΦ]
neC = [count(Z -> classify(Z, i) === :eC, Zs) for i in 1:nΦ]

# --------------------------------------------------------------------------
# 3. Validation against the hero
# --------------------------------------------------------------------------
hero = load(HERO_PATH)["res"]
Ic_mm_all = vec(getindex(findmax(hero.Js; dims = 2), 1))
Φ_mm = flux_of_B.(Ref(junctions["fitted_C3"].wireL), collect(hero.params.Brng))
keep = findall(i -> isfinite(Ic_mm_all[i]) && first(Φrng) <= Φ_mm[i] <= last(Φrng), eachindex(Φ_mm))
xs, ys = Φ_mm[keep], Ic_mm_all[keep]
p = sortperm(xs)
Ic_hero = [lininterp(xs[p], ys[p], Φ) for Φ in Φrng]

# Compare shapes: the hero is at 460 mK and this at T = 0, so an overall scale
# difference is expected. A shape difference is not.
#
# The comparison is reported twice. Near the lobe edges the gap closes, and a
# 460 mK calculation is suppressed there while a T = 0 one is not — so the edges
# disagree *for a physical reason* and a whole-window figure is dominated by it.
# The interior number is the one that tests the conventions; the edge number is
# reported alongside so the split is visible rather than hidden by a window
# choice.
sc = maximum(filter(isfinite, Ic)) / maximum(filter(isfinite, Ic_hero))
peak = maximum(filter(isfinite, Ic))
dev  = abs.(Ic .- sc .* Ic_hero) ./ peak

const INTERIOR = (0.70, 1.40)
inner = findall(Φ -> INTERIOR[1] <= Φ <= INTERIOR[2], Φrng)
edge  = setdiff(eachindex(Φrng), inner)
worst       = maximum(filter(isfinite, dev))
worst_inner = maximum(filter(isfinite, dev[inner]))
worst_edge  = isempty(edge) ? NaN : maximum(filter(isfinite, dev[edge]))

# The zero-energy crossings show up as steps in Ic. If the two models agree, the
# steps land at the same flux — a π error in the Zeeman term would move them.
steps(y) = [i for i in 2:length(y) if isfinite(y[i]) && isfinite(y[i-1]) &&
                                      abs(y[i] - y[i-1]) > 0.02 * peak]
steps_1D   = Φrng[steps(Ic)]
steps_hero = Φrng[steps(sc .* Ic_hero)]

@info "Σ_Z Ic vs hero" scale_factor=sc worst_all=worst worst_interior=worst_inner worst_edge=worst_edge

# --------------------------------------------------------------------------
# 4. Report + plot
# --------------------------------------------------------------------------
eC_frac = maximum(eC_peak) / maximum(hC_peak)
println("\n" * "="^72)
println("§3 DIAGNOSTIC — 1L lobe, Φ ∈ [$(round(first(Φrng);digits=3)), $(round(last(Φrng);digits=3))]")
println("="^72)
println("sectors swept                : ", Zs)
println("max Σ hC peak weight         : ", maximum(hC_peak))
println("max Σ eC peak weight         : ", maximum(eC_peak))
println("eC / hC peak ratio           : ", round(eC_frac; digits = 4))
println("eC sector count, range       : ", extrema(neC))
println("hC sector count, range       : ", extrema(nhC))
println("most negative signed eC      : ", minimum(eC_sgn))
println()
println("VALIDATION vs hero")
println("  scale factor               : ", round(sc; digits = 4))
println("  worst dev, interior $(INTERIOR): ", round(100worst_inner; digits = 2), " % of peak")
println("  worst dev, lobe edges       : ", round(100worst_edge; digits = 2), " % of peak")
println("  worst dev, whole window     : ", round(100worst; digits = 2), " % of peak")
println("  Ic steps, 1D  at Φ =        : ", round.(steps_1D; digits = 4))
println("  Ic steps, hero at Φ =       : ", round.(steps_hero; digits = 4))
println()

crossings = [i for i in 2:nΦ if neC[i] != neC[i-1]]
println("SECTOR CONVERSIONS (eC → hC) at Φ = ", round.(Φrng[crossings]; digits = 4))
print("  eC count by region          : ")
regions = vcat(1, crossings, nΦ + 1)
for k in 1:length(regions)-1
    i = regions[k]
    print(neC[i], k < length(regions)-1 ? " → " : "")
end
println()
println()
if eC_frac > 0.5
    println(">>> CASE: COMPETITION. eC weight rivals hC somewhere in the lobe.")
    println("    Pick three fluxes: hC-dominated, near crossover, eC-grown.")
elseif eC_frac < 0.05
    println(">>> CASE: SUBSTITUTION, eC a hairline.")
    println("    Pick two fluxes bracketing one crossing, plus one away.")
else
    println(">>> NEITHER CASE AS WRITTEN. eC/hC peak ratio is $(round(eC_frac; digits=3)) —")
    println("    a persistent minority, not a hairline and not comparable. The mechanism")
    println("    is substitution by discrete steps: eC sectors convert to hC with rising")
    println("    flux, and the total steps UP at each conversion. Build doc §1 and the §9")
    println("    caption describe the opposite direction; see the report.")
end
println("="^72)

mkpath(OUTDIR)
fig = Figure(size = (1100, 900))

ax1 = Axis(fig[1, 1]; ylabel = "Ic (a.u.)", title = "Σ_Z Ic vs hero (hero rescaled)")
lines!(ax1, Φrng, Ic; color = :black, linewidth = 2.5, label = "Σ_Z Ic (1D, T=0)")
lines!(ax1, Φrng, sc .* Ic_hero; color = :orchid, linestyle = :dash, linewidth = 2.5,
       label = "hero (mm, 460 mK) × $(round(sc; digits=3))")
axislegend(ax1; position = :lb)

ax2 = Axis(fig[2, 1]; ylabel = "peak weight (a.u.)", title = "peak weight by family")
lines!(ax2, Φrng, hC_peak; color = TURQ,  linewidth = 3, label = "Σ hC")
lines!(ax2, Φrng, eC_peak; color = CORAL, linewidth = 3, label = "Σ eC")
axislegend(ax2; position = :lt)

ax3 = Axis(fig[3, 1]; ylabel = "signed contribution", title = "contribution at the total's φ*")
hlines!(ax3, [0]; color = (:gray, 0.6))
lines!(ax3, Φrng, hC_sgn; color = TURQ,  linewidth = 3, label = "Σ hC at φ*")
lines!(ax3, Φrng, eC_sgn; color = CORAL, linewidth = 3, label = "Σ eC at φ*")
axislegend(ax3; position = :lb)

ax4 = Axis(fig[4, 1]; xlabel = "Φ/Φ₀", ylabel = "# sectors", title = "sector counts")
stairs!(ax4, Φrng, Float64.(nhC); color = TURQ,  linewidth = 3, label = "hC")
stairs!(ax4, Φrng, Float64.(neC); color = CORAL, linewidth = 3, label = "eC")
vlines!(ax4, Φrng[crossings]; color = (:black, 0.35), linestyle = :dot)
axislegend(ax4; position = :lt)

for a in (ax1, ax2, ax3, ax4); xlims!(a, first(Φrng), last(Φrng)); end
save(joinpath(OUTDIR, "cpr_diagnostic.png"), fig; px_per_unit = 2)
@info "wrote $(joinpath(OUTDIR, "cpr_diagnostic.png"))"
fig
