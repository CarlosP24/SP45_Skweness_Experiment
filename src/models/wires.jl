params_LDOS = Calc_Params(;
    Brng = range(0, 0.4, 101),
    ωrng = range(-0.3, 0, 101) .+ 1e-3im,
)

wires = Dict(
    "base_wire" => Params_mm(;
        R = 60,
        d = 8,
        a0 = 1,
        Δ0 = 0.2,
        α = 0,
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
    Δw::Float64 = 0
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

function add_wire(name::String)
    # Expected name shapes:
    #  - "C_left_mu=1.0_w=20_dw=3_tau=2.5_g=10"
    #  - "D_left_mu=5_w=10" (some params may be absent)
    parts = split(name, '_')
    if length(parts) < 2
        error("invalid wire name: $name")
    end

    base_key = string(parts[1], "_", parts[2])
    if !haskey(wires, base_key)
        error("unknown base wire: $base_key")
    end

    # collect k=v tokens from the remainder
    tokens = parts[3:end]
    params_map = Dict{String,String}()
    for t in tokens
        if occursin('=', t)
            k, v = split(t, '=', limit=2)
            params_map[k] = v
        end
    end

    # helper to parse numeric values robustly (accepts integers and floats)
    parse_num(s::AbstractString) = try
        parse(Float64, replace(s, ',' => '.'))
    catch
        error("could not parse numeric value '$s' in wire name '$name'")
    end

    # Build keyword args for Params_mm. The original code uses unicode fields
    # like `μ` and `τΓ`, so we construct a Dict{Symbol,Any} and splat it.
    kwargs = Dict{Symbol,Any}()
    if haskey(params_map, "mu")
        kwargs[Symbol("μ")] = parse_num(params_map["mu"])
    end
    if haskey(params_map, "w")
        kwargs[:w] = parse_num(params_map["w"])
    end
    if haskey(params_map, "tau")
        kwargs[Symbol("τΓ")] = parse_num(params_map["tau"])
    end
    if haskey(params_map, "g")
        kwargs[:g] = parse_num(params_map["g"])
    end
    # other possible params can be added here (e.g. num_mJ)

    # Create new Params_mm based on the base template and provided kwargs
    base = wires[base_key]
    newwire = isempty(kwargs) ? Params_mm(base) : Params_mm(base; kwargs...)

    # Δw (dw in the string) is stored on the Wire_System struct
    Δw = haskey(params_map, "dw") ? parse_num(params_map["dw"]) : 0.0

    ws = Wire_System(
        wire = newwire,
        params = params_LDOS,
        Δw = Δw,
    )

    # register and return
    wire_systems[name] = ws
end

function read_wire(name::String)
    parts = split(name, '_')
    # If the name looks like a wire descriptor (e.g. "C_left_..." or "D_right_...")
    # delegate to add_wire which will parse and register it.
    if length(parts) >= 2 && (parts[2] == "left" || parts[2] == "right")
        add_wire(name)
    end
    # Otherwise expect it to be a direct key in wire_systems
    if !haskey(wire_systems, name)
        @info "$name is not a wire"
    end
end

# μs = range(0, 50, step = 0.5)
# ws = range(5, 30, step = 1)
# Δws = range(0, 30, step = 1)
# τΓs = range(1, 5, step = 0.5)
# gs = range(0, 20, step = 5)
# nmJs = collect(round.(range(5, 10, length=length(μs)))) .|> Int
# pts = Iterators.product(zip(μs, nmJs), ws, Δws, τΓs, gs)

# wires_C_left = Dict(
#     ["C_left_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Wire_System(
#         wire = Params_mm(wires["C_left"]; μ = μ, w = w, num_mJ = nmJ, τΓ = τΓ, g = g),
#         params = params_LDOS,
#         Δw =  Δw
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# wires_C_right = Dict(
#     ["C_right_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Wire_System(
#         wire = Params_mm(wires["C_right"]; μ = μ, w = w, num_mJ = nmJ, τΓ = τΓ, g = g),
#         params = params_LDOS,
#         Δw = Δw
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# wires_D_left = Dict(
#     ["D_left_mu=$(μ)_w=$(w)" => Wire_System(
#         wire = Params_mm(wires["D_left"]; μ = μ, w = w, num_mJ = nmJ, τΓ = τΓ, g = g),
#         params = params_LDOS,
#         Δw = Δw
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# wires_D_right = Dict(
#     ["D_right_mu=$(μ)_w=$(w)" => Wire_System(
#         wire = Params_mm(wires["D_right"]; μ = μ, w = w, num_mJ = nmJ),
#         params = params_LDOS,
#         Δw = Δw
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# merge!(wire_systems, wires_C_left, wires_C_right, wires_D_left, wires_D_right)