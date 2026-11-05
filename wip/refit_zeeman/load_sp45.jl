# Load the SP45 source on all workers, exactly as src/main.jl does, plus the refit helpers.
# Local:  julia --project=. wip/refit_zeeman/<script>.jl   (SP45_WORKERS workers, default 4)
# Slurm:  the same command inside an allocation adds the Slurm tasks as workers.
using Distributed
if haskey(ENV, "SLURM_JOB_ID") && nprocs() == 1
    using SlurmClusterManager
    @time addprocs(SlurmManager(); exeflags = "--project=$(pwd())")
elseif nprocs() == 1
    addprocs(parse(Int, get(ENV, "SP45_WORKERS", "4")); exeflags = "--project=$(pwd())")
end

@everywhere begin
    using Pkg
    using Quantica, FullShell, JLD2, ProgressMeter, Parameters, Interpolations, Roots, QuadGK
    const SRC = joinpath(pwd(), "src")
    include(joinpath(SRC, "utilities.jl"))
    include(joinpath(SRC, "models", "params.jl"))
    include(joinpath(SRC, "models", "wires.jl"))
    include(joinpath(SRC, "models", "junctions.jl"))
    include(joinpath(SRC, "models", "cpr.jl"))
    include_all("builders")
    include_all("operators")
    include_all("parallelizers")
    include_all("calculations")
    include(joinpath(pwd(), "wip", "refit_zeeman", "ic_single_mode.jl"))
end
