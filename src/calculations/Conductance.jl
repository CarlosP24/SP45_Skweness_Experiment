function G0_finiteT(G::Function, kBT::Real; preω = 10, imω)
    df_dω(ω) = -1/(4*kBT) * sech(ω/(2*kBT))^2
    ωmax = preω * kBT
    result, _ = quadgk(ω -> G(ω) * df_dω(ω), -ωmax + imω*1im, -ωmax/2 + 10im*imω, 0 + imω*1im, ωmax/2 + 10im*imω, ωmax + imω*1im; rtol = 1e-6)
    return result
end
function calc_Conductance(name::String)
    junction = junctions[name]
    @unpack wireL, wireR, params, TN, kBT = junction
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

        # Build integration path
    bw = maximum([bandwidth(params_left), bandwidth(params_right)])
    itipL = get_itip(params_left)            
    itipR = get_itip(params_right)
    itip(x) = minimum([itipL(x), itipR(x)])
    ΩL = get_Ω(params_left)
    ΩR = get_Ω(params_right)
    Ω(x) = maximum([ΩL(x; θ = params_left.θ), ΩR(x; θ = params_right.θ)])
    ipath = Paths.polygon((mu, kBT; B = 0, _...) -> (-bw, -Ω(B),  -Ω(B)/2 + itip(B)*1im, 0) .+ imshift*1im)


    Gf = densitymatrix(g[attach_link["semi"]], (-bw, bw); omegamap = ω -> (; ω))

    function Gparalell(B)
        try
            if Ω(B) < 1e-4
                return NaN
            end
            return Gf(kBT; B, τ,)
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