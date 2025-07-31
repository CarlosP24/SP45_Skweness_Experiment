function filter_πs(φrng)
    φrng1 = φrng[findall(x -> x > 0 && x < π, φrng)]
    φrng2 = φrng[findall(x -> x > π && x < 2π, φrng)]
    deleteat!(φrng1, findall(x -> isapprox(x, π), φrng1))
    deleteat!(φrng2, findall(x -> isapprox(x, 2π), φrng2))
    return vcat(φrng1, φrng2)
end

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


function calc_Josephson(name::String)
    junction = junctions[name]
    @unpack wireL, wireR, params, j_params, TN = junction
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

    # Build Hamiltonian and Green's function
    hSM_left, hSC_left, params_left = build_cyl_mm(wireL )
    hSM_right, hSC_right, params_right = build_cyl_mm(wireR)
    g_right, g_left, g = greens_dict["semi"](hSC_left, hSC_right, params_left, params_right;)

    # Find τ from TN
    τrng = subdiv(0, 1, 100)
    Gτs = get_TN(hSM_left, hSM_right, params_left, params_right, "semi", τrng; B = 0)
    Gτs = Gτs ./ maximum(Gτs)
    Tτ = linear_interpolation(τrng, Gτs)            # gives TN as a function of τ
    τ = find_zeros(τ -> Tτ(τ) - TN, 0, 1) |> first
    println("Found τ = $τ for TN = $TN")

    # Build integration path
    bw = maximum([bandwidth(params_left), bandwidth(params_right)])
    itipL = get_itip(params_left)            
    itipR = get_itip(params_right)
    itip(x) = minimum([itipL(x), itipR(x)])
    ΩL = get_Ω(params_left)
    ΩR = get_Ω(params_right)
    Ω(x) = maximum([ΩL(x; θ = params_left.θ), ΩR(x; θ = params_right.θ)])
    ipath = Paths.polygon((mu, kBT; B = 0, _...) -> (-bw, -Ω(B),  -Ω(B)/2 + itip(B)*1im, 0) .+ imshift*1im)

    # Josephson integrator
    J = josephson(g[attach_link["semi"]], ipath; omegamap = ω -> (; ω), phases = φrng, atol, maxevals, order)

    function Jparalell(B)
        try
            return J(0; B, τ)
        catch
            return fill(NaN, length(φrng))
        end
    end

    # Calculate Josephson current
    Js = pfunction(
        Jparalell,
        [Brng]
    )

    Js = mapreduce(permutedims, vcat, Js)

    return Results(
        params = params2,
        junction = junction,
        Js = Js,
        path = path
    )
end