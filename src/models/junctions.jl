@with_kw struct J_Params
    imshift = 1e-4
    atol = 1e-7
    maxevals = 1e5
    order = 21
end

@with_kw struct Junction
    wireL::Params_mm = Params_mm()
    wireR::Params_mm = Params_mm()
    params::Calc_Params = Calc_Params()
    j_params::J_Params = J_Params()
    TN = 1
end


params_jos = Calc_Params(;
    Brng = range(0, 0.6, 101),
    φrng = range(0, 2π, 101),
)

#TNs = round.(10 .^ range(-4, 0, 5), sigdigits = 1)
TN = 1.0
μs = range(0, 5, step = 0.5)
ws = range(5, 30, step = 5)
nmJs = collect(round.(range(5, 10, length=length(μs)))) .|> Int
pts = Iterators.product(zip(μs, nmJs), ws)

junctions_C = Dict(
    ["C_mu=$(μ)_w=$(w)" => Junction(
        wireL = Params_mm(wires["C_left"]; μ = μ, w = w, num_mJ = nmJ),
        wireR = Params_mm(wires["C_right"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_jos,
        TN = TN,
    ) for ((μ, nmJ), w) in pts]
)

junctions_D = Dict(
    ["D_mu=$(μ)_w=$(w)" => Junction(
        wireL = Params_mm(wires["D_left"]; μ = μ, w = w, num_mJ = nmJ),
        wireR = Params_mm(wires["D_right"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_jos,
        TN = TN,
    ) for ((μ, nmJ), w) in pts]
)

junctions_dict = Dict(
    "junctions_C" => junctions_C,
    "junctions_D" => junctions_D
)

junctions = merge(junctions_C, junctions_D)