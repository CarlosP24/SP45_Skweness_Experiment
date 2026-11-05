# Step 3 of the devices-B/C fit: grid scan of I_c(B) with the single-sector builder and the physical
# Zeeman (src/calculations/Ic_single_mode.jl), α = 0, physical g between 0 and 20 (Carlos, 2026-11-05).
# Lead geometry from refit_BC/devices.jl. Free parameters as in the paper: μ, Γ_NS = τΓ Δ0, R_av = R − w/2
# (here w), T_N (an overall scale of I_c fitted in analyze.jl), plus g. τ is calibrated per (μ, w, τΓ) to
# T_N0 = 0.24 on the sector-summed normal conductance at B = 0 (g plays no role there). The kernel
# temperature is the T_e fitted per device in target.jl (B 604 → 600 mK, C 292 → 290 mK).
#
#   julia --project=. refit_BC/scan.jl <grid> [part] [tag]
# `part` = "r1,r2,.../K" keeps the candidates i with mod1(i, K) ∈ {r1, r2, ...} (interleaved, so every part
# samples the whole grid); output refit_BC/out/<grid>_<tag>.jld2. Merge the parts with analyze.jl.
include(joinpath(pwd(), "refit_BC", "load.jl"))

const NAME = isempty(ARGS) ? "smokeB" : ARGS[1]
const PART = length(ARGS) >= 2 ? ARGS[2] : "1/1"
const TAG = length(ARGS) >= 3 ? ARGS[3] : "all"

const GVALS = [0.0, 5.0, 10.0, 15.0, 20.0]
const MUS = [10.0, 15.0, 20.0, 25.0, 30.0]
const TAUS = [2.0, 2.5, 3.0, 3.5]
grids = Dict(
    # device B: R = 60.5 (lattice radius 60) → R_av = 52.5 … 42.5 nm; fields of the digitized map
    "B" => (dev = "B", T = 600.0, g = GVALS, w = [15.0, 20.0, 25.0, 30.0, 35.0], μ = MUS, τΓ = TAUS,
            Bs = collect(0.004:0.004:0.24), Zs = -9:9),
    # device C: R = 71 → R_av = 63.5 … 51 nm; lobes 0–3 (edge 3.5 at 0.414 T)
    "C" => (dev = "C", T = 290.0, g = GVALS, w = [15.0, 20.0, 25.0, 30.0, 35.0, 40.0], μ = MUS, τΓ = TAUS,
            Bs = collect(0.006:0.006:0.42), Zs = -12:12),
    # device C, extension past the edges the coarse optimum sits on (τΓ = 3.5, μ = 10–15; C_fixed_ranking.txt):
    # stronger coupling, fewer occupied sectors. |Z| ≤ 10 is ample at μ ≤ 10 meV (smokeC: −9…8 at μ = 30)
    "C2" => (dev = "C", T = 290.0, g = [0.0, 10.0, 20.0], w = [25.0, 30.0, 35.0, 40.0], μ = [2.5, 5.0, 7.5, 10.0],
             τΓ = [4.0, 5.0, 6.0], Bs = collect(0.006:0.006:0.42), Zs = -10:10),
    # device C, lead coherence length: device A's right lead (ξ_d = 92 nm) leaves only ≈ 0.3 Δ0 of shell gap
    # at C's lobe edges (R_LP = 75 nm), and every C candidate starts 1L too low. C has no Joule data, so ξ_d
    # (both leads) is scanned around the C optimum at a realistic g
    "C3" => (dev = "C", T = 290.0, g = [10.0], w = [25.0, 30.0, 35.0], μ = [10.0, 15.0, 20.0],
             τΓ = [2.5, 3.0, 3.5], ξd = [30.0, 45.0, 60.0], Bs = collect(0.006:0.006:0.42), Zs = -12:12),
    # device B, refinement around the best coarse candidate (g 10, w 30, μ 15, τΓ 3; B_fixed_ranking.txt)
    "B2" => (dev = "B", T = 600.0, g = [7.5, 10.0, 12.5], w = [27.5, 30.0, 32.5], μ = [12.5, 15.0, 17.5],
             τΓ = [2.75, 3.0, 3.25], Bs = collect(0.004:0.004:0.24), Zs = -9:9),
    # end-to-end pipeline test (scan → analyze → final → figure): one candidate on the full B grid
    "testB" => (dev = "B", T = 600.0, g = [15.0], w = [25.0], μ = [25.0], τΓ = [3.0],
                Bs = collect(0.004:0.004:0.24), Zs = -9:9),
    # pipeline and sector-range tests: the largest μ and R_av of each grid, extra sectors, few fields
    "smokeB" => (dev = "B", T = 600.0, g = [20.0], w = [15.0], μ = [30.0], τΓ = [3.0],
                 Bs = [0.02, 0.12, 0.2, 0.23], Zs = -12:12),
    "smokeC" => (dev = "C", T = 290.0, g = [20.0], w = [15.0], μ = [30.0], τΓ = [3.0],
                 Bs = [0.03, 0.12, 0.24, 0.35, 0.41], Zs = -15:15),
)
grid = grids[NAME]
cands = haskey(grid, :ξd) ?
    vec([(; μ, w, τΓ, g, ξd) for g in grid.g, w in grid.w, μ in grid.μ, τΓ in grid.τΓ, ξd in grid.ξd]) :
    vec([(; μ, w, τΓ, g) for g in grid.g, w in grid.w, μ in grid.μ, τΓ in grid.τΓ])
let (rs, K) = split(PART, '/')
    keep = parse.(Int, split(rs, ','))
    global cands = [c for (i, c) in enumerate(cands) if mod1(i, parse(Int, K)) in keep]
end

const DEV = grid.dev
const Bs = grid.Bs
const Zs = grid.Zs
const φrng = range(0, 2π, 101)
const TN0 = 0.24
const CHUNK = 4
@info "scan $NAME ($PART)" device = DEV ncands = length(cands) nB = length(Bs) nZ = length(Zs) T = grid.T workers = nworkers()

# τ for every distinct (μ, w, τΓ): one serial calibration per worker
keysτ = unique([(c.μ, c.w, c.τΓ) for c in cands])
τvals = pmap(keysτ) do (μ, w, τΓ)
    wL, wR = device_wires(DEV; μ, w, τΓ, g = 0.0)
    try
        find_τ_sm(wL, wR, TN0; Zs)
    catch e
        @warn "τ calibration failed" μ w τΓ exception = e
        NaN
    end
end
τs = Dict(zip(keysτ, τvals))
@info "τ calibrated" τs

chunks = collect(Iterators.partition(eachindex(Bs), CHUNK))
tasks = [(ic, ich) for ic in eachindex(cands) for ich in eachindex(chunks)]
t = @elapsed out = @showprogress pmap(tasks) do (ic, ich)
    c = cands[ic]
    τ = τs[(c.μ, c.w, c.τΓ)]
    wL, wR = device_wires(DEV; cand_kw(c)...)
    isfinite(τ) || return nothing
    js_reduced(wL, wR, τ, Bs[chunks[ich]]; T = grid.T, Zs, φrng)
end
@info "scan done" hours = t / 3600

nφ = length(filter_πs(φrng))
Jtot = fill(NaN, length(cands), length(Bs), nφ)
JZc = fill(NaN, length(cands), length(Bs), length(Zs))
IcZ = fill(NaN, length(cands), length(Bs), length(Zs))
for ((ic, ich), o) in zip(tasks, out)
    o === nothing && continue
    Jtot[ic, chunks[ich], :] .= o.Jtot
    JZc[ic, chunks[ich], :] .= o.JZc
    IcZ[ic, chunks[ich], :] .= o.IcZ
end
Ic = dropdims(maximum(Jtot; dims = 3); dims = 3)            # cand × B

mkpath(joinpath("refit_BC", "out"))
path = joinpath("refit_BC", "out", "$(NAME)_$(TAG).jld2")
jldsave(path; cands, τs = [τs[(c.μ, c.w, c.τΓ)] for c in cands], dev = DEV, Bs, Zs = collect(Zs),
    φs = filter_πs(φrng), Ic, Jtot, JZc, IcZ, TN0, T = grid.T, geometry = GEOMETRY[DEV],
    wall_hours = t / 3600, nworkers = nworkers())
@info "saved" path
