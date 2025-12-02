using Pkg, JLD2
include("utilities.jl")
ensure_pkgs(["Quantica", "FullShell", "Parameters", "ProgressMeter", "Interpolations", "Roots"])
@everywhere begin
    using Quantica, FullShell, JLD2
    using ProgressMeter, Parameters
    using Interpolations, Roots
    
    # Load code
    include("utilities.jl")
    include("models/params.jl")
    include("models/wires.jl")
    include("models/junctions.jl")
    include_all("builders")
    include_all("operators")
    include_all("parallelizers")
    include_all("calculations")
end

## Run 
input = ARGS[1]

read_wire(input)
read_junction(input)

if input in keys(wires_fit)
    @info "Fitting LDOS for wire: $input"
    res = fit_LDOS(input)
    @info "LDOS fitting completed for wire: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
    exit(0)
end

if input in keys(wire_systems)
    @info "Calculating LDOS for wire: $input"
    res = calc_LDOS(input)
    @info "LDOS calculation completed for wire: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
    exit(0)
end

if input in keys(junctions)
    @info "Calculating Josephson current for: $input"
    res = calc_Josephson(input)
    @info "Josephson current calculation completed for: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
    exit(0)
end

if endswith(input, "_cond")
    truekey = replace(input, "_cond" => "")
    if truekey in keys(junctions)
        @info "Calculating Conductance for: $truekey"
        res = calc_Conductance(truekey)
        @info "Conductance calculation completed for: $truekey"
        @info "Saving results to $(res.path)"
        save(res.path, "res", res)
        exit(0)
    end
end

if input in keys(junctions_fit)
    @info "Fitting Josephson current for: $input"
    res = fit_Josephson(input)
    @info "Josephson current fitting completed for: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
    exit(0)
end

@error "Key '$input' not found."