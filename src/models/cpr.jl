# m_J-resolved current-phase relations.
#
# ---------------------------------------------------------------------------
# WHY THERE IS A SECOND MODEL HERE AT ALL
#
# The hero figure's junction is built with `build_cyl_mm`: every m_J channel
# lives in a single Hamiltonian as a transverse lattice coordinate, and
# Quantica's `josephson` traces over the whole contact. So the per-m_J
# decomposition the CPR figure is *about* cannot be recovered from that object —
# `josephson` accepts a contact slice and nothing finer.
#
# `build_cyl` instead carries m_J as the model parameter `Z`, one sector per
# run. The decomposition is then exact by construction, which is how SP43
# produced the same panels.
#
# The price is that the two builders do not share conventions. Both differences
# below are silent — the code runs either way and returns a plausible-looking
# spectrum — so they are converted here, in one place, and nowhere else.
#
#   1. FLUX. `Params_mm` is parameterised by the axial field B (T); `Params` by
#      the normalised flux Φ/Φ₀. `build_cyl` applies cos(θ) itself, so Φ is
#      handed over without it.
#
#   2. ZEEMAN — the one that bites.
#
#          build_cyl_mm :  0.5 g μBΦ0 πoΦ0 B cos θ
#          build_cyl    :  0.5 g μBΦ0 Φ cos θ / (π RLP²)
#
#      Substituting Φ = πoΦ0 B RLP² into the second gives exactly the first
#      divided by π. The single-mode form is the physical one — it reduces to
#      0.5 g μB B — and the multimode form is π times larger.
#
#      g = 17 was *fitted* in the multimode model, so it is an effective value in
#      that convention, not a bare g-factor. Handing it straight to `build_cyl`
#      would shrink the Zeeman splitting by π and move every CdGM analog, which
#      is the entire subject of the figure. Hence g_1D = π g_mm below.
#
#      Do not "fix" this by changing g in the fit. The hero depends on the
#      multimode value, and this figure has to match the hero, not correct it.
#
# `verify_against_hero` (calculations/CPR.jl) is what actually demonstrates the
# two models agree. Treat it as part of the calculation, not as a test.
# ---------------------------------------------------------------------------

"""
    flux_of_B(p, B)

Normalised flux Φ/Φ₀ threading the Little–Parks area of wire `p` at axial field
`B` (T). Uses the wire's own `πoΦ0`, so the multimode value stays authoritative
(`Params` and `Params_mm` disagree in the 6th digit — 7 ppm, irrelevant to the
physics but not worth introducing gratuitously).
"""
flux_of_B(p, B) = p.πoΦ0 * B * (p.R + p.d / 2)^2

"Inverse of [`flux_of_B`](@ref)."
B_of_flux(p, Φ) = Φ / (p.πoΦ0 * (p.R + p.d / 2)^2)

"""
    single_mode(p::Params_mm; kw...)

The single-mode `Params` mirror of a fitted multimode wire, with the flux and
Zeeman conventions converted as described in this file's header.

`ishollow = true` reproduces the multimode geometry, which is a hollow core at a
single radius `Rav = R - w/2`; for the fitted wires that is 60 − 25/2 = 47.5 nm,
the radius the poster quotes for the CdGM analogs. The dome profile is switched
off (`Vmax = Vmin = 0`) because the multimode builder has none.
"""
function single_mode(p::Params_mm; kw...)
    return Params(;
        a0 = p.a0,
        R = p.R,
        w = p.w,
        d = p.d,
        Δ0 = p.Δ0,
        ξd = p.ξd,
        θ = p.θ,
        μ = p.μ,
        τΓ = p.τΓ,
        α = p.α,
        g = π * p.g,          # convention 2 — see header. Not a typo.
        shell = p.shell,
        ishollow = true,
        Vmax = 0, Vmin = 0,   # no dome profile in the multimode model
        L = 0,                # semi-infinite leads, as in the hero
        kw...
    )
end

"""
A junction whose current-phase relation is resolved by m_J sector.

`Zs` is the set of conserved angular-momentum sectors to sweep. It must cover
the multimode model's channel count, which is `|r[2]/a0| <= num_mJ`, i.e.
`-num_mJ:num_mJ` — 13 sectors for the fitted `num_mJ = 6`. Truncating it below
that silently drops current.

`source` names the multimode junction this mirrors, so the provenance sidecar and
`verify_against_hero` can find it.
"""
@with_kw struct CPR_System
    wireL::Params = Params()
    wireR::Params = Params()
    Zs = -6:6
    Φrng = subdiv(0.5, 1.5, 81)     # the 1L lobe — §3 of the build doc
    φrng = range(0, 2π, 101)        # identical to the hero's
    j_params::J_Params = J_Params()
    TN::Float64 = 0.24
    kB::Float64 = 8.61e-5           # meV/mK, matching the Junction struct
    T::Float64 = 0.0                # mK — deliberately 0, see below
    source::String = ""
    outdir::String = "data/theory/cpr"
end

cpr_systems = Dict{String, CPR_System}()

"""
    add_cpr(name, source; kw...)

Register an m_J-resolved CPR system mirroring the fitted multimode junction
`source`. Every wire parameter is *derived* from that junction rather than
retyped, so the two cannot drift apart — the "same parameters as the hero" claim
in the caption is then structural rather than a promise.

The one deliberate deviation is temperature. The hero runs at Te = 460 mK; this
runs at `T = 0`, because at 460 mK the steps from zero-energy crossings round off
and the decomposition stops being readable. `TN` stays at the fitted 0.24.
"""
function add_cpr(name::String, source::String; T = 0.0, kw...)
    haskey(junctions, source) || error("unknown source junction: $source")
    j = junctions[source]
    j.Δw == 0 || @warn "source junction $source has Δw = $(j.Δw); the single-mode " *
                       "mirror bakes w in at build time and cannot average over it"
    sys = CPR_System(;
        wireL = single_mode(j.wireL),
        wireR = single_mode(j.wireR),
        # num_mJ is a Float64 field on Params_mm (@deftype Float64); Z is an
        # angular-momentum index and is used as a Dict key, so keep it integral.
        Zs = -Int(j.wireL.num_mJ):Int(j.wireL.num_mJ),
        j_params = j.j_params,
        TN = j.TN,
        kB = j.kB,
        T = T,
        source = source,
        kw...
    )
    cpr_systems[name] = sys
    return sys
end

function read_cpr(name::String)
    startswith(name, "cpr_") || return
    source = name[5:end]
    haskey(junctions, source) && add_cpr(name, source)
    return
end

## Registered CPR systems
add_cpr("cpr_fitted_C3", "fitted_C3")
