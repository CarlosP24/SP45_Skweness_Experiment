using Pkg, JLD2
include("utilities.jl")
ensure_pkgs(["Quantica", "FullShell", "Parameters", "ProgressMeter"])
@everywhere begin
    using Quantica, FullShell
    using ProgressMeter, Parameters
    
    # Load code
    include("utilities.jl")
    include_all("builders")
    include_all("models")
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
else
    @error "Key '$input' not found."
end
