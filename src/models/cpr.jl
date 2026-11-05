# m_J-resolved current-phase relations.
#
# ---------------------------------------------------------------------------
# All SP45 wires are single-sector `build_cyl` models (src/models/wires.jl): m_J
# is the model parameter `Z`, and the junction's current is Σ_Z of the sector
# currents. `calc_Josephson` sums the sectors right away; this file keeps them
# apart, which is what the CPR figure is about.
#
# History. The published curves were first computed with the multimode builder
# `build_cyl_mm`, where every m_J channel is a transverse lattice coordinate and
# the decomposition cannot be recovered. That builder also writes the Zeeman
# term π times too large (0.5 g μBΦ0 πoΦ0 B instead of 0.5 g μB B), so its
# "g = 17" is a physical g = 17π ≈ 53.4. The single-sector wires now carry the
# physical g directly (`g_fit` in wires.jl), and no conversion is needed here.
#
# `verify_against_hero` (calculations/CPR.jl) checks that Σ_Z of this
# calculation reproduces the junction's I_c.
# ---------------------------------------------------------------------------

"""
    flux_of_B(p, B)

Normalised flux Φ/Φ₀ threading the Little–Parks area of wire `p` at axial field
`B` (T), with the wire's own `πoΦ0` (the multimode value, see wires.jl).
"""
flux_of_B(p, B) = p.πoΦ0 * B * (p.R + p.d / 2)^2

"Inverse of [`flux_of_B`](@ref)."
B_of_flux(p, Φ) = Φ / (p.πoΦ0 * (p.R + p.d / 2)^2)

"""
    single_mode(p::Params; kw...)

The wire of the CPR calculation: the junction's own single-sector wire, with
optional overrides. Kept as a function so the CPR system is derived from the
junction rather than retyped. g is the physical g, as everywhere.
"""
single_mode(p::Params; kw...) = Params(p; kw...)

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
        # Z is an angular-momentum index and is used as a Dict key: the wire's own range
        Zs = j.wireL.Zs,
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
