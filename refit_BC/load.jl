# Load the SP45 source on all workers, exactly as src/main.jl does, plus the device definitions.
# Local:  julia --project=. refit_BC/<script>.jl   (SP45_WORKERS workers, default 4)
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
    include_all("calculations")             # includes Ic_single_mode.jl
    include(joinpath(pwd(), "refit_BC", "devices.jl"))

    """
    single-sector left and right wires of device `dev` for one fit candidate; `ξd` (nm), when given,
    replaces the coherence length of both leads (devices without Joule characterisation)
    """
    function device_wires(dev; μ, w, τΓ, g, α = 0.0, ξd = nothing)
        G = GEOMETRY[dev]
        ξL, ξR = ξd === nothing ? (G.ξdL, G.ξdR) : (ξd, ξd)
        wL = sm_wire(; R = G.R, d = G.d, w, ξd = ξL, θ = deg2rad(G.θL), Δ0 = G.Δ0, μ, τΓ, g, α)
        wR = sm_wire(; R = G.R, d = G.d, w, ξd = ξR, θ = deg2rad(G.θR), Δ0 = G.Δ0, μ, τΓ, g, α)
        return wL, wR
    end
    "the device_wires keywords of a candidate (its ξd only if the grid varies it)"
    cand_kw(c) = haskey(c, :ξd) ? (; μ = c.μ, w = c.w, τΓ = c.τΓ, g = c.g, ξd = c.ξd) :
                                  (; μ = c.μ, w = c.w, τΓ = c.τΓ, g = c.g)
end
