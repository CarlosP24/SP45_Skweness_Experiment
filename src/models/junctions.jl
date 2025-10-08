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

params_jos = Calc_Params(;
    Brng = range(0, 0.4, 201),
    φrng = range(0, 2π, 101),
)


# Registry for dynamically created junctions
junctions = Dict{String, Junction}()


function add_junction(name::String)
    """
    Parse a junction name and create/register a Junction.

    Expected name shapes (tokens after base are optional):
      - "C_mu=1.0_w=20_dw=3_tau=2.5_g=10"
      - "D_mu=5_w=10"

    The first token is the base: "C" or "D". This function uses
    wires["<base>_left"] and wires["<base>_right"] as templates
    for the left/right Params_mm and applies any provided numeric overrides.
    """

    parts = split(name, '_')
    if length(parts) < 1
        error("invalid junction name: $name")
    end

    base = parts[1]
    if !(base in ("C", "D"))
        error("unknown junction base: $base (expected 'C' or 'D')")
    end

    # collect k=v tokens from the remainder
    tokens = parts[2:end]
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
        error("could not parse numeric value '$s' in junction name '$name'")
    end

    # Build keyword args for Params_mm. Map mu->μ, tau->τΓ, w->w, g->g
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

    # Create wire Params_mm for left and right from templates
    left_key = string(base, "_left")
    right_key = string(base, "_right")
    if !haskey(wires, left_key) || !haskey(wires, right_key)
        error("missing wire templates for junction base: $base")
    end

    baseL = wires[left_key]
    baseR = wires[right_key]

    wireL = isempty(kwargs) ? Params_mm(baseL) : Params_mm(baseL; kwargs...)
    wireR = isempty(kwargs) ? Params_mm(baseR) : Params_mm(baseR; kwargs...)

    # Δw (dw in the string) is stored on the Junction struct
    Δw = haskey(params_map, "dw") ? parse_num(params_map["dw"]) : 0.0

    j = Junction(
        wireL = wireL,
        wireR = wireR,
        params = params_jos,
        Δw = Δw,
        j_params = J_Params(),
        TN = 0.1,
    )

    # register and return
    junctions[name] = j
end

function read_junction(name::String)
    parts = split(name, '_')
    # If the name looks like a junction descriptor (e.g. "C_..." or "D_...")
    if length(parts) >= 1 && (parts[1] == "C" || parts[1] == "D")
        add_junction(name)
    end
    # Otherwise, not a recognized junction descriptor
    @info "$name is not a valid junction descriptor"
end




#########################################################################





# fitC = try 
#     load("data/Js_fit/base_C.jld2")["res"].fitted_params
# catch e 
#     load("../data/Js_fit/base_C.jld2")["res"].fitted_params
# end

# μs = range(0, 50, step = 0.5)
# ws = range(5, 30, step = 1)
# Δws = range(5, 30, step = 1)
# τΓs = range(1, 5, step = 0.5)
# gs = range(0, 20, step = 5)
# nmJs = collect(round.(range(5, 10, length=length(μs)))) .|> Int

# pts = Iterators.product(zip(μs, nmJs), ws, Δws, τΓs, gs)


# junctions_C = Dict(
#     ["C_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Junction(
#         wireL = Params_mm(wires["C_left"]; μ, w, τΓ, g),
#         wireR = Params_mm(wires["C_right"]; μ, w, τΓ, g),
#         params = params_jos,
#         Δw = Δw,
#         TN = 0.1,
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# junctions_D = Dict(
#     ["D_mu=$(μ)_w=$(w)_dw=$(Δw)_tau=$(τΓ)_g=$(g)" => Junction(
#         wireL = Params_mm(wires["D_left"]; μ, w, τΓ, g),
#         wireR = Params_mm(wires["D_right"]; μ, w, τΓ, g),
#         params = params_jos,
#         Δw = Δw,
#         TN = 0.1,
#     ) for ((μ, nmJ), w, Δw, τΓ, g) in pts]
# )

# junctions_dict = Dict(
#     "junctions_C" => junctions_C,
#     "junctions_D" => junctions_D
# )

# junctions = merge(junctions_C, junctions_D)