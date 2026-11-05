# Josephson current and LDOS from the m_J-independent (single-sector) builder `build_cyl`, summed over
# the conserved sectors Z, with the physical Zeeman energy ½ g μ_B B cos θ.
#
# Why: the multimode builder `build_cyl_mm` (used for every published SP45 curve) writes the Zeeman term
# π times too large (wip/refit_zeeman/REPORT.md). The junction is cylindrically symmetric and both leads
# have equal R and d, hence the same fluxoid number, so J_z is conserved through the device: the multimode
# Hamiltonian is block diagonal in m_J and Σ_Z of the sector currents is the same current. At the
# published Zeeman (g_1D = 17π) the sector sum reproduces data/Js/fitted_C3.jld2 within 0.5 % of the
# peak (wip/refit_zeeman/validate.txt; refit_BC/validate.jl repeats the check with these functions).
#
# Parameters are the FullShell `Params` of a hollow core (tubular-core approximation of the paper):
# radius R_av = R − w/2, no dome profile, semi-infinite leads, a0 = 1 nm, α = 0 unless given. The flux
# Φ/Φ0 = π R_LP² B / Φ0 (R_LP = R + d/2) is the model parameter; build_cyl applies cos θ itself.
# Everything here is serial inside a worker, so it can be called from a pmap over candidates.

"""
    sm_wire(; R, d, w, ξd, θ, Δ0, μ, τΓ, g, α = 0, kw...)

Single-sector wire with the physical g-factor `g` (θ in radians).
"""
sm_wire(; R, d, w, ξd, θ, Δ0, μ, τΓ, g, α = 0.0, kw...) =
    Params(; a0 = 1, R, d, w, ξd, θ, Δ0, μ, τΓ, g, α, shell = "Ballistic", ishollow = true,
             Vmax = 0, Vmin = 0, L = 0, kw...)

"Φ/Φ0 through the shell of `p` at axial field B (T), with the multimode convention πoΦ0 of the published fit"
flux_sm(p::Params, B) = 1.519257e-3 * B * (p.R + p.d / 2)^2

"""
    build_J_sm(wireL, wireR; j_params, φrng)

Josephson integrator of the single-sector junction (as calc_CPR) and the smaller of the two shell gaps
as a function of flux. `J(kBT; Φ, Z, τ)` returns J_S at the phases `φs` (0, π, 2π removed).
"""
function build_J_sm(wireL::Params, wireR::Params; j_params = J_Params(), φrng = range(0, 2π, 101))
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
    js_reduced(wireL, wireR, τ, Bs; T, Zs, φrng, j_params, kB)

Current of the junction at the fields `Bs`, reduced to what the fit needs, so a worker returns kilobytes:
  Jtot[iB, iφ] = Σ_Z J_Z(φ)               (the total current-phase relation)
  JZc[iB, iZ]  = J_Z(φ*), φ* = argmax_φ Jtot (each sector's share of I_c)
  IcZ[iB, iZ]  = max_φ |J_Z(φ)|            (each sector's own critical current)
NaN where either shell is gapless (Ω < 1 µeV) or the integration fails. T in mK.
"""
function js_reduced(wireL::Params, wireR::Params, τ, Bs; T, Zs, φrng = range(0, 2π, 101),
        j_params = J_Params(), kB = 8.61e-5)
    (; J, Ωmin, φs) = build_J_sm(wireL, wireR; j_params, φrng)
    nB, nZ, nφ = length(Bs), length(Zs), length(φs)
    Jtot = fill(NaN, nB, nφ); JZc = fill(NaN, nB, nZ); IcZ = fill(NaN, nB, nZ)
    for (i, B) in enumerate(Bs)
        Φ = flux_sm(wireL, B)
        Ωmin(Φ) < 1e-3 && continue
        JZ = fill(NaN, nZ, nφ)
        for (k, Z) in enumerate(Zs)
            try
                JZ[k, :] .= J(kB * T; Φ, Z, τ)
            catch e
                @warn "J failed at B = $B, Z = $Z" exception = e
            end
        end
        tot = vec(sum(JZ; dims = 1))
        Jtot[i, :] .= tot
        IcZ[i, :] .= vec(maximum(abs.(JZ); dims = 2))
        all(isfinite, tot) && (JZc[i, :] .= JZ[:, argmax(tot)])
    end
    return (; Jtot, JZc, IcZ, φs)
end

"""
    find_τ_sm(wireL, wireR, TN; Zs, τrng)

Coupling τ at which the normal-state conductance Σ_Z G_Z(τ) at B = 0 (shell replaced by the broadening
Γ = τΓ Δ0, as `get_TN`) reaches the fraction `TN` of its maximum over τ: the multimode `find_τ` of
calculations/Josephson.jl, sector by sector. Serial.
"""
function find_τ_sm(wireL::Params, wireR::Params, TN; Zs, τrng = subdiv(0, 1, 100))
    hL, _, pL = build_cyl(wireL)
    hR, _, pR = build_cyl(wireR)
    _, _, g = greens_dict["semi"](add_Δ0(hL, pL), add_Δ0(hR, pR), pL, pR;)
    G = conductance(g[attach_link["semi"], attach_link["semi"]])
    Gτ = [sum(G(0; τ, Φ = 0, Z) for Z in Zs) for τ in τrng]
    Tτ = linear_interpolation(τrng, Gτ ./ maximum(Gτ))
    return find_zeros(τ -> Tτ(τ) - TN, 0, 1) |> first
end

"""
    ldos_sm(wire, Bs, ωs; Zs)

LDOS at the end cell of the semi-infinite lead `wire` (as calc_LDOS), summed over orbitals and over the
sectors `Zs`: a length(Bs) × length(ωs) matrix. Serial.
"""
function ldos_sm(wire::Params, Bs, ωs; Zs)
    _, hSC, p = build_cyl(wire)
    g_right = hSC |> greenfunction(GS.Schur(boundary = 0))
    ρ = ldos(g_right[cells = (-1,)])
    out = zeros(length(Bs), length(ωs))
    for (i, B) in enumerate(Bs), (j, ω) in enumerate(ωs)
        Φ = flux_sm(p, B)
        acc = 0.0
        for Z in Zs
            acc += try
                sum(ρ(ω; Φ, ω, Z))
            catch
                NaN
            end
        end
        out[i, j] = acc
    end
    return out
end
