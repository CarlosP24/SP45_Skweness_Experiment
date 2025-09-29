params_LDOS = Calc_Params(;
    Brng = range(0, 0.4, 201),
    ωrng = range(-0.3, 0, 101) .+ 1e-4im,
)

wires = Dict(
    "base_wire" => Params_mm(;
        R = 60,
        d = 8,
        a0 = 2,
        Δ0 = 0.2,
        α = 10,
        g = 20,
        # Parameters to fit 
        w = 20,
        τΓ = 3,
        μ = 0,
        num_mJ = 4,
        shell = "Ballistic",
    )
)


wires["C_left"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 55,
    #θ = 8 * π / 180,
    θ = 0.1,
    Δ0 = 0.204,
)

wires["C_right"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 92,
    #θ = 5 * π / 180,
    θ = 0.1,
    Δ0 = 0.204
)

wires["D_left"] = Params_mm(wires["base_wire"];
    R = 52,
    ξd = 68,
    θ = 0,
    Δ0 = 0.203
)

wires["D_right"] = Params_mm(wires["base_wire"];
    R = 51,
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

μs = range(0, 5, step = 0.5)
ws = range(5, 30, step = 5)
nmJs = collect(round.(range(5, 10, length=length(μs)))) .|> Int
pts = Iterators.product(zip(μs, nmJs), ws)

wires_C_left = Dict(
    ["C_left_mu=$(μ)_w=$(w)" => Wire_System(
        wire = Params_mm(wires["C_left"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_LDOS
    ) for ((μ, nmJ), w) in pts]
)

wires_C_right = Dict(
    ["C_right_mu=$(μ)_w=$(w)" => Wire_System(
        wire = Params_mm(wires["C_right"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_LDOS
    ) for ((μ, nmJ), w) in pts]
) 

wires_D_left = Dict(
    ["D_left_mu=$(μ)_w=$(w)" => Wire_System(
        wire = Params_mm(wires["D_left"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_LDOS
    ) for ((μ, nmJ), w) in pts]
)

wires_D_right = Dict(
    ["D_right_mu=$(μ)_w=$(w)" => Wire_System(
        wire = Params_mm(wires["D_right"]; μ = μ, w = w, num_mJ = nmJ),
        params = params_LDOS
    ) for ((μ, nmJ), w) in pts]
)

merge!(wire_systems, wires_C_left, wires_C_right, wires_D_left, wires_D_right)