params_LDOS = Calc_Params(;
    Brng = range(0, 0.4, 401),
    ωrng = range(-0.3, 0, 401) .+ 1e-3im,
    #ωrng = [0] .+ 1e-3im
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
        num_mJ = 6,
        shell = "Ballistic",
    )
)


wires["C_left"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 55,
    θ = 8 * π / 180,
    Δ0 = 0.204,
)

wires["C_right"] = Params_mm(wires["base_wire"];
    R = 60,
    ξd = 92,
    θ = 5 * π / 180,
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

## Fitted wires 

function add_fit(name::String, key)
    add_wire("C_left_"*name)
    add_wire("C_right_"*name)
    wire_systems["fitted_"*key*"_left"] = wire_systems["C_left_"*name]
    wire_systems["fitted_"*key*"_right"] = wire_systems["C_right_"*name]
    delete!(wire_systems, "C_left_"*name)
    delete!(wire_systems, "C_right_"*name)
end

C = "mu=25_w=20_dw=0_tau=2.5_g=17"
C3 = "mu=20_w=25_tau=2.5_g=17_T=550"
add_fit(C, "C")
add_fit(C3, "C3")