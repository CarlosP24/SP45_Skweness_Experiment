function calc_Conductance(name::String)
    junction = junctions[name]
    @unpack wireL, wireR, params, TN = junction
    @unpack Brng, ωrng, outdir = params

    imω = imag(ωrng[1])

    Brng = collect(Brng)
    deleteat!(Brng, findall(x -> isapprox(x, 0), Brng))
    params2 = Calc_Params(params; Brng)


    # Setup output path
    path = "$(outdir)/Conductance/$(name).jld2"
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

    Gf = conductance(g[attach_link["semi"], attach_link["semi"]]; nambu = true)

    function Gparalell(B)
        try
            return Gf(0.0 .+ imω*1im; ω = 0.0 .+ imω*1im, τ, B)
        catch e 
            @warn "Error calculating G at B=$B: $e"
            return NaN
        end 
    end

    Gs = pfunction(Gparalell, [Brng])

    return Results(
        params = params2,
        junction = junction,
        Gs = Gs,
        path = path,
    )
end