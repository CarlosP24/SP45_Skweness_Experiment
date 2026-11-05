# Critical current of the SP45 junction summed over conserved m_J sectors, with the
# single-mode builder `build_cyl` (physical Zeeman ½gμ_B B), for the Zeeman-corrected refit
# (Carlos, 2026-11-05). The junction is cylindrically symmetric and both leads have the same
# R and d, hence the same fluxoid number, so J_z is conserved through the device and the
# multimode Hamiltonian of `build_cyl_mm` is block diagonal in m_J: Σ_Z J_Z(φ) is the same
# current, without the π-too-large Zeeman of `build_cyl_mm`. Validated against the published
# multimode fit at g_1D = 17π (validate.txt: ≤ 0.5 % of the peak).
#
# Load after the SP45 src (see load_sp45.jl). Everything here is new; nothing in src/ changes.

"""
    sm_wire(p::Params_mm; g, kw...)

Single-mode mirror of the multimode wire `p` with the *physical* g-factor `g`. The keyword
replaces the π·p.g that `single_mode` (src/models/cpr.jl) uses to reproduce the published fit.
"""
sm_wire(p::Params_mm; g = p.g, kw...) = single_mode(p; g = g, kw...)

"""
    build_J(wireL, wireR; j_params, φrng)

Josephson integrator of the single-mode junction (as calc_CPR) and the smaller of the two
shell gaps as a function of flux. `J(kBT; Φ, Z, τ)` returns J_S at the phases `φs`.
"""
function build_J(wireL::Params, wireR::Params; j_params = J_Params(), φrng = range(0, 2π, 101))
    @unpack imshift, atol, maxevals, order = j_params
    φs = filter_πs(φrng)
    _, hSC_left, params_left = build_cyl(wireL)
    _, hSC_right, params_right = build_cyl(wireR)
    _, _, g = greens_dict["semi"](hSC_left, hSC_right, params_left, params_right;)
    bw = maximum([bandwidth(params_left), bandwidth(params_right)])
    itipL, itipR = get_itip(params_left), get_itip(params_right)
    itip(x) = minimum([itipL(x), itipR(x)])
    ΩL, ΩR = get_Ω(params_left), get_Ω(params_right)
    Ω(x) = maximum([ΩL(x), ΩR(x)])
    Ωmin(x) = minimum([ΩL(x), ΩR(x)])
    ipath = Paths.polygon((mu, kBT; Φ = 0, _...) -> (-bw, -Ω(Φ), -Ω(Φ)/2 + itip(Φ)*1im, 0) .+ imshift*1im)
    J = josephson(g[attach_link["semi"]], ipath; omegamap = ω -> (; ω), phases = φs, atol, maxevals, order)
    return (; J, Ωmin, φs)
end

"""
    js_serial(wireL, wireR, τ, Φs; T, Zs, φrng, kB)

Serial version for use inside a worker: returns a (length(Φs), length(Zs), length(φs)) array
of J_S, NaN where either shell is gapless (Ω < 1 µeV) or the integration fails. T in mK.
"""
function js_serial(wireL::Params, wireR::Params, τ, Φs; T = 460.0, Zs = -6:6,
        φrng = range(0, 2π, 101), j_params = J_Params(), kB = 8.61e-5)
    (; J, Ωmin, φs) = build_J(wireL, wireR; j_params, φrng)
    out = fill(NaN, length(Φs), length(Zs), length(φs))
    for (i, Φ) in enumerate(Φs)
        Ωmin(Φ) < 1e-3 && continue
        for (k, Z) in enumerate(Zs)
            try
                out[i, k, :] .= J(kB * T; Φ, Z, τ)
            catch e
                @warn "J failed at Φ = $Φ, Z = $Z" exception = e
            end
        end
    end
    return out
end

"""
    ic_single_mode(wireL, wireR, τ, Φs; T, Zs, φrng)

Σ_Z J_Z(φ) at each flux in `Φs`, maximised over φ, parallel over (Φ, Z) with `pcpr`.
Returns (Ic, Js, φs) with Js a Dict(Z => length(Φs) × length(φs) matrix).
"""
function ic_single_mode(wireL::Params, wireR::Params, τ, Φs; T = 460.0, Zs = -6:6,
        φrng = range(0, 2π, 101), j_params = J_Params(), kB = 8.61e-5)
    (; J, Ωmin, φs) = build_J(wireL, wireR; j_params, φrng)
    function Jsector(Φ, Z)
        Ωmin(Φ) < 1e-3 && return fill(NaN, length(φs))
        return J(kB * T; Φ, Z, τ)
    end
    Js = pcpr(Jsector, Φs, Zs, length(φs))
    Jtot = sum(values(Js))
    Ic = vec(maximum(Jtot; dims = 2))
    return (; Ic, Js, φs)
end
