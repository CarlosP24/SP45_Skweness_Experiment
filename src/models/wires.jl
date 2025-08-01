params_LDOS = Calc_Params(;
    Brng = range(0, 0.6, 101),
    ωrng = range(-0.3, 0, 101) .+ 1e-4im,
)

wires = Dict(
    "base_wire" => Params_mm(;
        R = 60,
        d = 8,
        a0 = 2,
        Δ0 = 0.2,
        α = 0,
        g = 10,
        # Parameters to fit 
        w = 20,
        τΓ = 1,
        μ = 0,
        num_mJ = 5,
        shell = "Ballistic",
    )
)


wires["C_left"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 55,
    θ = 8 * π / 180,
    Δ0 = 0.204
)

wires["C_right"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 92,
    θ = 5 * π / 180,
    Δ0 = 0.204
)

wires["D_left"] = Params_mm(wires["base_wire"];
    R = 56,
    ξd = 68,
    θ = 0,
    Δ0 = 0.203
)

wires["D_right"] = Params_mm(wires["base_wire"];
    R = 55,
    ξd = 70,
    θ = 0,
    Δ0 = 0.2
)

@with_kw struct Wire_System
    wire::Params_mm = Params_mm()
    params::Calc_Params = Calc_Params()
end

wire_systems = Dict(
    [
        key => Wire_System(
            wire = value, 
            params = params_LDOS
        )
        for (key, value) in wires
    ]
)