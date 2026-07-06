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
    TN = 0.1
    kB::Float64 = 8.61e-5   # eV/K
    T::Float64 = 0.0        # mK
end

params_jos = Calc_Params(;
    #Brng = vcat(collect(0.01:0.005:0.02),collect(0.05:0.005:0.08), collect(0.18:0.005:0.22), collect(0.28:0.005:0.32)),
    Brng = range(0, 0.38, 201),
    φrng = range(0, 2π, 101),
)


# Registry for dynamically created junctions
junctions = Dict{String, Junction}()

## Fitted junctions 


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
        @info "$name is not a junction"
        return
    end

    base = parts[1]
    if !(base in ("C", "D"))
        @info "$name has unknown junction base: $base (expected 'C' or 'D')"
        return
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
        @info "missing wire templates for junction base: $base"
        return
    end

    baseL = wires[left_key]
    baseR = wires[right_key]

    wireL = isempty(kwargs) ? Params_mm(baseL) : Params_mm(baseL; kwargs...)
    wireR = isempty(kwargs) ? Params_mm(baseR) : Params_mm(baseR; kwargs...)

    # Δw (dw in the string) is stored on the Junction struct
    Δw = haskey(params_map, "dw") ? parse_num(params_map["dw"]) : 0.0
    T = haskey(params_map, "T") ? parse_num(params_map["T"]) : 0.0
    TN = haskey(params_map, "TN") ? parse_num(params_map["TN"]) : 0.24

    j = Junction(
        wireL = wireL,
        wireR = wireR,
        params = params_jos,
        Δw = Δw,
        j_params = J_Params(),
        TN = TN,
        T = T,
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


## Fitted junctions
add_junction("C_mu=25_w=20_dw=0_tau=2.5_g=17_T=550_TN=0.19")
add_junction("C_mu=20_w=29_tau=3.1_g=20_T=550")
add_junction("C_mu=20_w=25_tau=2.5_g=17_T=550")
junctions["fitted_C"] = junctions["C_mu=25_w=20_dw=0_tau=2.5_g=17_T=550_TN=0.19"]
junctions["fitted_C2"] = junctions["C_mu=20_w=29_tau=3.1_g=20_T=550"]
junctions["fitted_C3"] = junctions["C_mu=20_w=25_tau=2.5_g=17_T=550"]
delete!(junctions, "C_mu=25_w=20_dw=0_tau=2.5_g=17_T=550_TN=0.19")
delete!(junctions, "C_mu=20_w=29_tau=3.1_g=20_T=550")
delete!(junctions, "C_mu=20_w=25_tau=2.5_g=17_T=550")

#########################################################################