function fit_Josephson(name::String)
    junction = junctions_fit[name]
    @unpack wireL, wireR, params, j_params, TN = junction
    @unpack µrng, τΓrng, wrng, grng, φrng, Bs1, Bs2, Idrop1, skew2, outdir = params
    @unpack imshift, atol, maxevals, order = j_params

    # Delete pathological points
    φrng = filter_πs(φrng)
    params2 = Fit_Params(params;  φrng)

    if skew2 == "right"
        Bs2 = reverse(Bs2)
    end

    # Setup output path
    path = "$(outdir)/Js_fit/$(name).jld2"
    mkpath(dirname(path))

    # Build Hamiltonian and Green's function
    hSM_left, hSC_left, params_left = build_cyl_mm(wireL)
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

    #ipath = Paths.polygon((mu, kBT; B = 0, _...) -> (-bw, 0) .+ imshift*1im)
    # Josephson integrator
    J = josephson(g[attach_link["semi"]], ipath; omegamap = ω -> (; ω), phases = φrng, atol, maxevals, order)

    function Jparalell(µ, τΓ, w, ga)
        try
            JB(B) = J(0; B, µ, τΓ, w, g = ga, τ)
            # Criteria 1: Check I_c drop between lobes
            Ic01 = maximum.(JB.(Bs1))
            if Idrop1 * Ic01[1] >= Ic01[2]
                # Criteria 2: check skewness of 2nd lobe
                Ic2 = maximum.(JB.(Bs2))
                if Ic2[1] > Ic2[2]
                    return true
                end
            end
            return false
        catch e
            @warn "Error in Jparalell at µ=$(µ), τΓ=$(τΓ), w=$(w), g=$(g)" exception = e
            return false
        end
    end

    # Calculate Josephson current
    Jfits = pfunction(
        Jparalell,
        [µrng, τΓrng, wrng, grng]
    )

    pts = collect(Iterators.product(µrng, τΓrng, wrng, grng))
    true_indices = findall(Jfits)
    fitted_params = pts[true_indices]

    return Results(
        params = params2,
        junction = junction,
        fitted_params = fitted_params,
        path = path
    )
end