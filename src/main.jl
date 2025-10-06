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

if input in keys(wire_systems)
    @info "Calculating LDOS for wire: $input"
    res = calc_LDOS(input)
    @info "LDOS calculation completed for wire: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
elseif input in keys(junctions)
    @info "Calculating Josephson current for: $input"
    res = calc_Josephson(input)
    @info "Josephson current calculation completed for: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
elseif input in keys(junctions_fit)
    @info "Fitting Josephson current for: $input"
    res = fit_Josephson(input)
    @info "Josephson current fitting completed for: $input"
    @info "Saving results to $(res.path)"
    save(res.path, "res", res)
else
    @error "Key '$input' not found."
end
