params_LDOS = Calc_Params(;
    Brng = range(0, 0.4, 401),
    ωrng = range(-0.3, 0, 401) .+ 1e-3im,
    #ωrng = [0] .+ 1e-3im
)

# Wires are single m_J sectors of FullShell's `build_cyl` (hollow core at R_av = R − w/2, no dome
# profile, semi-infinite leads). The junction conserves J_z (both leads have the same R and d), so every
# observable is the sum over the sectors `Zs`; Zs = -num_mJ:num_mJ is the truncation of the former
# multimode model (`build_cyl_mm`, |m_J| ≤ num_mJ), which these wires replace.
#
# g is the PHYSICAL Landé factor: `build_cyl` writes the Zeeman energy as ½ g μ_B B cos θ. The
# multimode builder wrote it π times larger, so the published fit's "g = 17" is a physical
# g = 17π ≈ 53.4. See `g_fit` below and README.md.
#
# πoΦ0 is set to the multimode value (it differs from the `Params` default in the 6th digit), so the
# flux Φ/Φ0 = πoΦ0 · B · R_LP² is exactly that of the published calculations.
wires = Dict(
    "base_wire" => Params(;
        R = 60,
        d = 8,
        a0 = 1,
        Δ0 = 0.2,
        α = 0,
        g = 20,
        # Parameters to fit
        w = 20,
        τΓ = 3,
        μ = 0,
        num_mJ = 6,
        Zs = -6:6,
        shell = "Ballistic",
        ishollow = true,
        Vmax = 0, Vmin = 0,
        L = 0,
        πoΦ0 = 1.519257e-3,
    )
)


wires["C_left"] = Params(wires["base_wire"];
    R = 60,
    ξd = 55,
    θ = 8 * π / 180,
    Δ0 = 0.204,
)

wires["C_right"] = Params(wires["base_wire"];
    R = 60,
    ξd = 92,
    θ = 5 * π / 180,
    Δ0 = 0.204
)

wires["D_left"] = Params(wires["base_wire"];
    R = 52,
    ξd = 68,
    θ = 0,
    Δ0 = 0.203
)

wires["D_right"] = Params(wires["base_wire"];
    R = 51,
    ξd = 70,
    θ = 0,
    Δ0 = 0.2
)

@with_kw struct Wire_System
    wire::Params = Params()
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

    # Build keyword args for Params. The original code uses unicode fields
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

    # Create new Params based on the base template and provided kwargs
    base = wires[base_key]
    newwire = isempty(kwargs) ? Params(base) : Params(base; kwargs...)

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

# Physical g of the fits. The published fit was done with the multimode builder at "g = 17", i.e. a
# physical g = 17π = 53.41; the paper quotes it as the fitted g ≈ 53, and `g_fit` is that value.
# (fitted_C, an earlier fit not in the paper, keeps its physical 17π.)
const g_fit = parse(Float64, get(ENV, "SP45_GFIT", "53"))
C = "mu=25_w=20_dw=0_tau=2.5_g=$(17π)"
C3 = "mu=20_w=25_tau=2.5_g=$(g_fit)_T=550"
add_fit(C, "C")
add_fit(C3, "C3")