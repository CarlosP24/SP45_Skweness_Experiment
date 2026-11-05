# Zeeman-corrected refit of the Paya:26c junction, step 2: grid scan of I_c(B).
# Single-mode sectors with the physical Zeeman (ic_single_mode.jl), physical g ∈ [15, 20]
# (Carlos, 2026-11-05). Everything else as the published fit `fitted_C3`: R = 60, d = 8 nm,
# ξ_L = 55, ξ_R = 92 nm, tilts 8° and 5°, Δ0 = 0.204 meV, Ballistic shell, a0 = 1 nm, α = 0,
# T = 460 mK. τ is calibrated per (μ, w, τΓ) to T_N = 0.24 with the hero's multimode normal-state
# conductance at B = 0 (Zeeman-free, so unaffected by the bug); T_N is refitted afterwards as an
# overall scale of I_c (D11 of ch5-plan: a transparency rescaling).
#
#   julia --project=. wip/refit_zeeman/scan.jl <name>        (name selects a grid below)
# Output: wip/refit_zeeman/out/<name>.jld2
include(joinpath(pwd(), "wip", "refit_zeeman", "load_sp45.jl"))

const NAME = isempty(ARGS) ? "scan1" : ARGS[1]

grids = Dict(
    # first coarse scan
    "scan1" => (g = [15.0, 17.5, 20.0], w = [15.0, 20.0, 25.0, 30.0, 35.0],
                μ = [10.0, 15.0, 20.0, 25.0, 30.0], τΓ = [2.0, 2.5, 3.0, 3.5]),
    # g above 20 (Carlos, 2026-11-05: find the minimal g with a left skew in lobe 2), around the
    # region where lobe 1 fits; 17π reproduces the published (multimode) Zeeman as an anchor
    "scan2" => (g = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 17π], w = [25.0, 30.0, 35.0],
                μ = [10.0, 15.0, 20.0, 25.0], τΓ = [2.5, 3.0, 3.5]),
    # where the strong lobe-2 left skew switches on between g = 45 (−0.004 T) and 17π (−0.025 T),
    # at the fitted radius and the published/best lobe-1 couplings
    "gsweep" => (g = [44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0], w = [25.0],
                 μ = [20.0, 25.0], τΓ = [2.5, 3.0]),
    # Rashba SOC at the physical g (Carlos, 2026-11-05): does α ≈ 20–30 meV nm give the lobe-2 left
    # skew that the α = 0 fit needed g ≈ 53 for? Both signs of α relative to the code's g.
    "alpha" => (g = [15.0, 17.5, 20.0], w = [25.0, 30.0], μ = [15.0, 20.0, 25.0], τΓ = [2.5, 3.0],
                α = [-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0]),
    # pipeline test: two candidates, few fields
    "smoke" => (g = [17.0], w = [25.0], μ = [20.0], τΓ = [2.5, 3.0]),
    "smoke_alpha" => (g = [17.0], w = [25.0], μ = [20.0], τΓ = [3.0], α = [0.0, 20.0]),
)
grid = grids[NAME]
αs = haskey(grid, :α) ? grid.α : [0.0]            # Rashba α (meV nm), constant at R_av; 0 if absent
cands = vec([(; μ, w, τΓ, g, α) for g in grid.g, w in grid.w, μ in grid.μ, τΓ in grid.τΓ, α in αs])

# Optional split across clusters: ARGS[2] = "r1,r2,.../K" keeps the candidates i with
# mod1(i, K) ∈ {r1, r2, ...} (interleaved, so every part samples the whole grid). Output tag
# ARGS[3] (default "all"). Merge the parts with analyze.jl.
const PART = length(ARGS) >= 2 ? ARGS[2] : "1/1"
const TAG = length(ARGS) >= 3 ? ARGS[3] : "all"
let (rs, K) = split(PART, '/')
    keep = parse.(Int, split(rs, ','))
    global cands = [c for (i, c) in enumerate(cands) if mod1(i, parse(Int, K)) in keep]
end

const Bs = startswith(NAME, "smoke") ? collect(range(0.05, 0.35, 4)) : collect(range(0.005, 0.38, 76))
const Zs = -9:9
const φrng = range(0, 2π, 101)
const TN0 = 0.24
const CHUNK = 4

hero = junctions["fitted_C3"]
@everywhere hero = junctions["fitted_C3"]
@info "scan $NAME" ncands = length(cands) nB = length(Bs) nZ = length(Zs) workers = nworkers()

# τ for every distinct (μ, w, τΓ, α): g plays no role at B = 0
keysτ = unique([(c.μ, c.w, c.τΓ, c.α) for c in cands])
τs = Dict{NTuple{4,Float64},Float64}()
for (μ, w, τΓ, α) in keysτ
    j = Junction(hero; wireL = Params_mm(hero.wireL; μ, w, τΓ, α), wireR = Params_mm(hero.wireR; μ, w, τΓ, α), TN = TN0)
    τs[(μ, w, τΓ, α)] = try
        find_τ(j)
    catch e
        @warn "τ calibration failed" μ w τΓ α exception = e
        NaN
    end
end
@info "τ calibrated" τs

@everywhere function scan_task(c, τ, Bchunk, Zs, φrng)
    pL = Params_mm(hero.wireL; μ = c.μ, w = c.w, τΓ = c.τΓ, α = c.α)
    pR = Params_mm(hero.wireR; μ = c.μ, w = c.w, τΓ = c.τΓ, α = c.α)
    wL, wR = sm_wire(pL; g = c.g), sm_wire(pR; g = c.g)
    Φs = flux_of_B.(Ref(pL), Bchunk)
    isfinite(τ) || return fill(NaN, length(Bchunk), length(Zs), length(filter_πs(φrng)))
    return js_serial(wL, wR, τ, Φs; T = hero.T, Zs, φrng, kB = hero.kB)
end

chunks = collect(Iterators.partition(eachindex(Bs), CHUNK))
tasks = [(ic, ich) for ic in eachindex(cands) for ich in eachindex(chunks)]
t = @elapsed out = @showprogress pmap(tasks) do (ic, ich)
    c = cands[ic]
    scan_task(c, τs[(c.μ, c.w, c.τΓ, c.α)], Bs[chunks[ich]], Zs, φrng)
end
@info "scan done" hours = t / 3600

nφ = length(filter_πs(φrng))
J = fill(NaN, length(cands), length(Bs), length(Zs), nφ)
for ((ic, ich), o) in zip(tasks, out)
    J[ic, chunks[ich], :, :] .= o
end
Jtot = dropdims(sum(J; dims = 3); dims = 3)                 # cand × B × φ
IcZ = dropdims(maximum(abs.(J); dims = 4); dims = 4)         # cand × B × Z, per-sector |J| max
Ic = dropdims(maximum(Jtot; dims = 3); dims = 3)             # cand × B

mkpath(joinpath("wip", "refit_zeeman", "out"))
jldsave(joinpath("wip", "refit_zeeman", "out", "$(NAME)_$(TAG).jld2");
    cands, τs = [τs[(c.μ, c.w, c.τΓ, c.α)] for c in cands], Bs, Zs = collect(Zs), φs = filter_πs(φrng),
    Ic, Jtot, IcZ, TN0, T = hero.T, wall_hours = t / 3600, nworkers = nworkers())
@info "saved" path = joinpath("wip", "refit_zeeman", "out", "$(NAME)_$(TAG).jld2")
