# filter_πs moved to src/utilities.jl — calc_CPR and the diagnostic need it too.

function add_Δ0(h, params)
    @unpack τΓ, Δ0 = params
    ΣS! = @onsite!((o, r; ) -> o - Δ0 * τΓ * im * σ0τ0)
    return h |> ΣS!
end

function get_TN(hleft, hright, params_left, params_right, gs, τrng; kw...)
    hc_left = add_Δ0(hleft, params_left)
    hc_right = add_Δ0(hright, params_right)
    g_right, g_left, g = greens_dict[gs](hc_left, hc_right, params_left, params_right;)
    G = conductance(g[attach_link[gs], attach_link[gs]])
    Gτ = pmap(τrng) do τ
        G(0; τ, kw...)
    end
    return reshape(Gτ, size(τrng)...)
end


"""
    find_τ(junction)

Coupling τ at which the normal-state conductance at B = 0, summed over the m_J sectors of the
junction, reaches the fraction `TN` of its maximum over τ (the shell enters only as the broadening
Γ = τΓ Δ0, see `get_TN`). There is no Zeeman term at B = 0, so this is the τ of the multimode
calibration of the published fit (0.86241 for fitted_C3).
"""
function find_τ(junction::Junction)
    @unpack wireL, wireR, TN = junction
    hSM_left, _, params_left = build_cyl(wireL)
    hSM_right, _, params_right = build_cyl(wireR)
    τrng = subdiv(0, 1, 100)
    Gτs = sum(get_TN(hSM_left, hSM_right, params_left, params_right, "semi", τrng; Φ = 0, Z) for Z in wireL.Zs)
    Gτs = Gτs ./ maximum(Gτs)
    Tτ = linear_interpolation(τrng, Gτs)            # gives TN as a function of τ
    return find_zeros(τ -> Tτ(τ) - TN, 0, 1) |> first
end

"""
    calc_Josephson(name)

Josephson current J_S(B, φ) of junction `name`: the sum over the conserved m_J sectors `Zs` of the
single-sector junctions (`build_cyl`), at the fields `Brng` (flux Φ = πoΦ0 B R_LP², `build_cyl`
applies cos θ itself). The sum over sectors is taken before any maximum over φ.
"""
function calc_Josephson(name::String)
    junction = junctions[name]
    @unpack wireL, wireR, params, Δw, j_params, TN, kB, T, τ_from = junction
    @unpack Brng, φrng, outdir = params
    @unpack imshift, atol, maxevals, order = j_params

    # Delete pathological points
    Brng = collect(Brng)
    deleteat!(Brng, findall(x -> isapprox(x, 0), Brng))
    φrng = filter_πs(φrng)
    params2 = Calc_Params(params; Brng, φrng)

    # Setup output path
    path = "$(outdir)/Js/$(name).jld2"
    mkpath(dirname(path))

    # Find τ from TN (of this junction, or of `τ_from` if given)
    τ = isempty(τ_from) ? find_τ(junction) : find_τ(junctions[τ_from])
    println("Found τ = $τ for TN = $TN" * (isempty(τ_from) ? "" : " of $τ_from"))

    Φ = get_Φ(wireL)                         # B (T) -> Φ/Φ0 through the shell
    Zs = collect(wireL.Zs)

    # The core width is a lattice property of build_cyl (R_av = R - w/2), so each w is its own build
    ws = [wireL.w]
    if Δw != 0
        push!(ws, wireL.w + Δw)
    end

    Js = zeros(length(Brng), length(φrng))
    for w in ws
        wL, wR = Params(wireL; w), Params(wireR; w)

        # Build Hamiltonian and Green's function
        _, hSC_left, params_left = build_cyl(wL)
        _, hSC_right, params_right = build_cyl(wR)
        g_right, g_left, g = greens_dict["semi"](hSC_left, hSC_right, params_left, params_right;)

        # Build integration path (a function of the flux)
        bw = maximum([bandwidth(params_left), bandwidth(params_right)])
        itipL = get_itip(params_left)
        itipR = get_itip(params_right)
        itip(x) = minimum([itipL(x), itipR(x)])
        ΩL = get_Ω(params_left)
        ΩR = get_Ω(params_right)
        Ω(x) = maximum([ΩL(x), ΩR(x)])
        Ωmin(x) = minimum([ΩL(x), ΩR(x)])
        ipath = Paths.polygon((mu, kBT; Φ = 0, _...) -> (-bw, -Ω(Φ),  -Ω(Φ)/2 + itip(Φ)*1im, 0) .+ imshift*1im)

        # Josephson integrator
        J = josephson(g[attach_link["semi"]], ipath; omegamap = ω -> (; ω), phases = φrng, atol, maxevals, order)

        function Jparalell(B, Z)
            try
                if Ωmin(Φ(B)) < 1e-3
                    return fill(NaN, length(φrng))
                end
                return J(kB*T; Φ = Φ(B), Z, τ)
            catch e
                @warn "Error in Jparalell at B=$(B), Z=$(Z)" exception = e
                return fill(NaN, length(φrng))
            end
        end

        # Calculate Josephson current, sector by sector, and sum the sectors
        JBZ = pfunction(Jparalell, [Brng, Zs])
        Js .+= mapreduce(permutedims, vcat, dropdims(sum(JBZ, dims = 2), dims = 2))
    end

    return Results(
        params = params2,
        junction = junction,
        Js = Js,
        path = path
    )
end
