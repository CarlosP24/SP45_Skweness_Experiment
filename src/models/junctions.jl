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
    Δw::Float64 = 0
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
    Brng = range(0, 0.4, 201),
    φrng = range(0, 2π, 101),
)

# fitC = try 
#     load("data/Js_fit/base_C.jld2")["res"].fitted_params
# catch e 
#     load("../data/Js_fit/base_C.jld2")["res"].fitted_params
# end

μs = range(0, 5, step = 0.5)
ws = range(5, 30, step = 5)
Δws = range(5, 30, step = 5)
τΓs = range(1, 5, step = 1)
gs = range(0, 20, step = 5)
nmJs = collect(round.(range(5, 10, length=length(μs)))) .|> Int

pts = Iterators.product(zip(μs, nmJs), ws, Δws, τΓs, gs)


junctions_C = Dict(
    ["C_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Junction(
        wireL = Params_mm(wires["C_left"]; μ, w, τΓ, g),
        wireR = Params_mm(wires["C_right"]; μ, w, τΓ, g),
        params = params_jos,
        Δw = Δw,
        TN = 0.1,
    ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
)

junctions_D = Dict(
    ["D_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Junction(
        wireL = Params_mm(wires["D_left"]; μ, w, τΓ, g),
        wireR = Params_mm(wires["D_right"]; μ, w, τΓ, g),
        params = params_jos,
        Δw = Δw,
        TN = 0.1,
    ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
)

junctions_dict = Dict(
    "junctions_C" => junctions_C,
    "junctions_D" => junctions_D
)

junctions = merge(junctions_C, junctions_D)