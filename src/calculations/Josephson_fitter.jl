function fit_Josephson(name::String)
    junction = junctions_fit[name]
    @unpack wireL, wireR, params, j_params, TN = junction
    @unpack µrng, τΓrng, wrng, grng, φrng, Δws, Bs01, Bs1, Bs2, Idrop01, Idrop01l, Ibump1, skew2, outdir, prevfit = params
    @unpack imshift, atol, maxevals, order = j_params

    fitprev = false
    if typeof(prevfit) == String
        try
        fitprev = load(prevfit)["res"].fitted_params
        catch 
        end
    end

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

    function Jparalell(µ, τΓ, w, Δw, ga)
        try
            JB = if Δw == 0
                B -> J(0; B, µ, τΓ, w = w, g = ga, τ)
            else
                B -> sum([J(0; B, µ, τΓ, w = w + dw, g = ga, τ) for dw in [0, Δw]])
            end
            # Criteria 1: 2nd lobe is skewed
            Ic2 = maximum.(JB.(Bs2))
            if Ic2[1] < Ic2[2]
                return false
            end
            # Criteria 2: Check I_c bumps in lobe 1
            Ic1 = maximum.(JB.(Bs1))
            if Ibump1 * Ic1[1] >= Ic1[2]
                return false
            end
            # Criteria 3: Check I_c drop between lobes, left side
            Ic01 = maximum.(JB.(Bs01))
            if Idrop01l * Ic01[1] <= Ic01[2]
                return false
            end
            # Criteria 4: Check I_c drop between lobes, right side
            if Idrop01 * Ic01[1] <= Ic01[3]
                return false
            end
            return true
        catch e
            @warn "Error in Jparalell at µ=$(µ), τΓ=$(τΓ), w=$(w), g=$(g)" exception = e
            return false
        end
    end

    if fitprev != false
        println("Evaluating only previously fitted points")
        µallowed = getindex.(fitprev, 1)
        τΓallowed = getindex.(fitprev, 2)
        wallowed = getindex.(fitprev, 3)
        Δwallowed = getindex.(fitprev, 4)
        gallowed = getindex.(fitprev, 4)
        µrng = intersect(µrng, µallowed)
        τΓrng = intersect(τΓrng, τΓallowed)
        wrng = intersect(wrng, wallowed)
        Δws = intersect(Δws, Δwallowed)
        grng = intersect(grng, gallowed)
    end

    # Calculate Josephson current
    Jfits = pfunction(
        Jparalell,
        [µrng, τΓrng, wrng, Δws, grng]
    )

    pts = collect(Iterators.product(µrng, τΓrng, wrng, Δws, grng))
    true_indices = findall(Jfits)
    fitted_params = pts[true_indices]

    return Results(
        params = params2,
        junction = junction,
        fitted_params = fitted_params,
        path = path
    )
end