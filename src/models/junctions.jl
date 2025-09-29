@with_kw struct J_Params
    imshift = 1e-4
    atol = 1e-7
    maxevals = 1e4
    order = 21
end

@with_kw struct Junction
    wireL::Params_mm = Params_mm()
    wireR::Params_mm = Params_mm()
    params::Calc_Params = Calc_Params()
    j_params::J_Params = J_Params()
    TN = 1
end

@with_kw struct Junction_fit
    wireL::Params_mm = Params_mm()
    wireR::Params_mm = Params_mm()
    params::Fit_Params = Fit_Params()
    j_params::J_Params = J_Params()
    TN = 0.1
end

junctions_fit = Dict(
    ["base_C" => Junction_fit(
        wireL = Params_mm(wires["C_left"]),
        wireR = Params_mm(wires["C_right"]),
        params = Fit_Params(),
        j_params = J_Params(),
    )]
)





#########################################################################

params_jos = Calc_Params(;
    Brng = range(0, 0.4, 101),
    φrng = range(0, 2π, 101),
)

#TNs = round.(10 .^ range(-4, 0, 5), sigdigits = 1)
TN = 0.1
µrng = range(0, 50, step = 1)
τΓrng = range(1, 5, step = 0.5)
wrng = range(5, 60, step = 5)
grng = range(0, 20, step = 5)

pts = Iterators.product(μrng, τΓrng, wrng, grng)

junctions_C = Dict(
    ["C_mu=$(μ)_w=$(w)_tau=$(τΓ)_g=$(g)" => Junction(
        wireL = Params_mm(wires["C_left"]; μ = μ, w = w, τΓ = τΓ, g = g),
        wireR = Params_mm(wires["C_right"]; μ = μ, w = w, τΓ = τΓ, g = g),
        params = params_jos,
        TN = TN,
    ) for (µ, τΓ, w, g) in pts]
)

junctions_D = Dict(
    ["D_mu=$(μ)_w=$(w)_tau=$(τΓ)_g=$(g)" => Junction(
        wireL = Params_mm(wires["D_left"]; μ = μ, w = w, τΓ = τΓ, g = g),
        wireR = Params_mm(wires["D_right"]; μ = μ, w = w, τΓ = τΓ, g = g),
        params = params_jos,
        TN = TN,
    ) for (μ, τΓ, w, g) in pts]
)

junctions_dict = Dict(
    "junctions_C" => junctions_C,
    "junctions_D" => junctions_D
)

junctions = merge(junctions_C, junctions_D)