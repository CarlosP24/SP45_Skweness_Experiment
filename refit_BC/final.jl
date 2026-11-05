# Step 5 of the devices-B/C fit: the curves of the figure for the selected candidates — I_c(B) on a fine
# field grid and the LDOS of the left and right leads (end of a semi-infinite wire, as Fig. 3(a, b)),
# both summed over the m_J sectors with the physical Zeeman.
#   julia --project=. refit_BC/final.jl <grid name> [T_e mode] [rank]
# Takes candidate `rank` (default 1) of refit_BC/out/<grid>_<mode>_scores.jld2 (analyze.jl).
# Output: refit_BC/out/final_<grid>_<mode>_<rank>.jld2
include(joinpath(pwd(), "refit_BC", "load.jl"))

const NAME = isempty(ARGS) ? "B" : ARGS[1]
const MODE = length(ARGS) >= 2 ? ARGS[2] : "fixed"
const RANK = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

sc = load(joinpath("refit_BC", "out", "$(NAME)_$(MODE)_scores.jld2"))
i = sc["order"][RANK]
c = sc["cands"][i]
DEV, TK = sc["dev"], sc["TK"]
const Zs = DEV == "B" ? (-9:9) : (-12:12)
const Bmax = DEV == "B" ? 0.24 : 0.42
const DB = parse(Float64, get(ENV, "SP45_FINAL_DB", "0.002"))       # field step of I_c (T)
const Bfine = collect(DB:DB:Bmax)
const BLDOS = collect(range(0, Bmax, parse(Int, get(ENV, "SP45_FINAL_NLDOS", "241"))))
const ωs = collect(range(-0.3, 0, 201)) .+ 1e-3im        # as params_LDOS (401 points there)
@info "final curves" NAME MODE RANK cand = c cost = sc["costs"][i] T = TK workers = nworkers()

wL, wR = device_wires(DEV; cand_kw(c)...)
τ = find_τ_sm(wL, wR, 0.24; Zs)
@info "τ" τ

chunks = collect(Iterators.partition(Bfine, 2))
t1 = @elapsed res = pmap(ch -> js_reduced(wL, wR, τ, collect(ch); T = TK, Zs), chunks)
Jtot = reduce(vcat, [r.Jtot for r in res]); JZc = reduce(vcat, [r.JZc for r in res]); IcZ = reduce(vcat, [r.IcZ for r in res])
Ic = vec(maximum(Jtot; dims = 2))
@info "I_c done" minutes = t1 / 60

t2 = @elapsed begin
    LL = pmap(b -> ldos_sm(wL, [b], ωs; Zs), BLDOS)
    LR = pmap(b -> ldos_sm(wR, [b], ωs; Zs), BLDOS)
end
LDOS_L = reduce(vcat, LL); LDOS_R = reduce(vcat, LR)
@info "LDOS done" minutes = t2 / 60

out = joinpath("refit_BC", "out", "final_$(NAME)_$(MODE)_$(RANK).jld2")
jldsave(out; cand = c, dev = DEV, T = TK, τ, s = sc["s"][i], Te = sc["Te"][i], cost = sc["costs"][i],
    Bs = Bfine, Ic, Jtot, JZc, IcZ, Zs = collect(Zs), BLDOS, ωs, LDOS_L, LDOS_R, geometry = GEOMETRY[DEV])
@info "saved" out
