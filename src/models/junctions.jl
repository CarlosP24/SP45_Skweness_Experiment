@with_kw struct J_Params
    imshift = 1e-6im
    atol = 1e-7
    maxevals = 1e7
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

TNs = 10 .^ range(-4, 0, 5)

junctions_C = Dict(
    ["C_$(TN)" => Junction(
        wireL = wires["C_left"],
        wireR = wires["C_right"],
        params = params_jos,
        TN = TN,
    ) for TN in TNs]
)

junctions_D = Dict(
    ["D_$(TN)" => Junction(
        wireL = wires["D_left"],
        wireR = wires["D_right"],
        params = params_jos,
        TN = TN,
    ) for TN in TNs]
)

junctions_dict = Dict(
    "junctions_C" => junctions_C,
    "junctions_D" => junctions_D
)

junctions = merge(junctions_C, junctions_D)