# Step 4 of the devices-B/C fit: score every scanned candidate against the device's target.
#   julia --project=plots refit_BC/analyze.jl <grid name> [T_e mode: fixed | free] [target]
# For each candidate the overall scale s of I_c (a transparency rescaling, T_N ≈ 0.24 s) is fitted, and
#   cost = Σ_{lobes n ≥ 1} rms[R_th − R_0/R_N]_n + λ0L · rms[(s I_c − I_sw)/I_sw(0)]_{0L}
# with R_th = I_0(I_c/0.0209 T_e/2)^(−2) (theo_exp_core.jl) and I_c in nA through the published 74.4.
# Every lobe n ≥ 1 has the same weight: no emphasis on, and no target for, the lobe-2 skewness.
# T_e mode "fixed" (default): the T_e fitted on the data (target.jl), as the paper does; the scan's kernel
# temperature is that value rounded. "free": T_e is fitted too (it only relates I_c to R_0).
# `target` (default: the scan's device) selects refit_BC/out/target_<target>.jld2, e.g. Cline.
# Skew diagnostics (not in the cost): centroid of each lobe's R_0/R_N dip, weight max(1 − R, 0), minus the
# lobe centre n·P; positive = skewed right.
# Outputs in refit_BC/out/: <name>_<mode>_ranking.txt, _per_g.txt, _scores.jld2; refit_BC/figures/<name>_<mode>_top.png
using JLD2, CairoMakie, Bessels, Statistics, Printf
include(joinpath(pwd(), "refit_BC", "devices.jl"))

const NAME = isempty(ARGS) ? "B" : ARGS[1]
const MODE = length(ARGS) >= 2 ? ARGS[2] : "fixed"
const TARGET = length(ARGS) >= 3 ? ARGS[3] : nothing
const CONV = 74.4          # code units → nA (plots/plotters/theo_exp_core.jl)
const λ0L = 0.5
const S_GRID = exp.(range(log(0.1), log(4.0), 300))
const TE_GRID = exp.(range(log(0.5), log(2.0), 61))    # T_e / T_kernel in mode "free"
const MARGIN = 0.002       # T, fields this close to a lobe edge are left out of the R cost

outdir = joinpath("refit_BC", "out")
files = sort([joinpath(outdir, f) for f in readdir(outdir) if startswith(f, NAME * "_") && endswith(f, ".jld2") &&
              !occursin("_fixed_", f) && !occursin("_free_", f)])
isempty(files) && error("no scan output for $NAME in $outdir")
parts = load.(files)
d = first(parts)
all(p -> p["Bs"] == d["Bs"], parts) || error("parts of $NAME use different field grids")
cands = reduce(vcat, [p["cands"] for p in parts])
Ic = reduce(vcat, [p["Ic"] for p in parts])
τs = reduce(vcat, [p["τs"] for p in parts])
Bs, DEV, TK = d["Bs"], d["dev"], d["T"]
println("merged $(length(files)) file(s), $(length(cands)) candidates of device $DEV, kernel T = $TK mK")

TGT = TARGET === nothing ? DEV : TARGET
tg = load(joinpath(outdir, "target_$TGT.jld2"))
const TE0 = tg["Te"]                                  # T_e of the target (mK)
edges = tg["edges"]
P = edges[2] - edges[1]                              # LP period (T)
centres = [(edges[n] + edges[n+1]) / 2 for n in 1:length(edges)-1]

"linear interpolation of y(x) at xq, NaN outside or across NaN"
function linterp(x, y, xq)
    i = searchsortedlast(x, xq)
    (i < 1 || i >= length(x)) && return NaN
    t = (xq - x[i]) / (x[i+1] - x[i])
    return (1 - t) * y[i] + t * y[i+1]
end
RAH(Ic_nA, Te) = besseli0(Ic_nA / (0.0209 * Te) / 2)^(-2)

# target fields: lobes n ≥ 1 inside the theory range, and the resolved 0L peak currents
Bmax = maximum(Bs)
# n for a field inside lobe n ≥ 1 bounded by two measured edges; nothing in 0L or past the last edge
function lobe_of(b)
    n = findlast(e -> b > e, edges)
    return (n === nothing || n >= length(edges)) ? nothing : n
end
mR = [b <= Bmax && lobe_of(b) !== nothing && minimum(abs.(b .- edges)) > MARGIN && isfinite(r)
      for (b, r) in zip(tg["B"], tg["r0"])]
Bt, r0, lob = tg["B"][mR], tg["r0"][mR], [lobe_of(b) for b in tg["B"][mR]]
lobes = sort(unique(lob))
m0 = [0.0 < b < edges[1] - 0.008 && isfinite(i) for (b, i) in zip(tg["B"], tg["Isw"])]
B0L, Isw0L = tg["B"][m0], tg["Isw"][m0]
Isw00 = maximum(Isw0L)
println("R cost: $(length(Bt)) fields in lobes $(lobes); 0L: $(length(B0L)) fields, I_sw(0) = $(round(Isw00; digits = 1)) nA")

"centroid of the R dip of lobe n (weight max(1 − R, 0)) minus the lobe centre; NaN without a dip"
function dip_offset(B, R, n)
    m = [lobe_of(b) == n && isfinite(r) for (b, r) in zip(B, R)]
    w = max.(1 .- R[m], 0)
    sum(w) < 1e-6 && return NaN
    return sum(B[m] .* w) / sum(w) - centres[n]
end
off_data = Dict(n => dip_offset(Bt, r0, n) for n in lobes)

function score(icurve)
    nan = (cost = Inf, s = NaN, Te = NaN, c0 = NaN, cL = Dict{Int,Float64}(), off = Dict{Int,Float64}())
    icB = [linterp(Bs, icurve, b) for b in Bt]
    ic0 = [linterp(Bs, icurve, b) for b in B0L]
    okB, ok0 = findall(isfinite, icB), findall(isfinite, ic0)
    (length(okB) < 0.9 * length(icB) || length(ok0) < 0.9 * length(ic0)) && return nan
    tes = MODE == "free" ? TE_GRID .* TE0 : [TE0]
    best = nan
    for s in S_GRID
        c0 = sqrt(mean(((CONV * s .* ic0[ok0]) .- Isw0L[ok0]) .^ 2)) / Isw00
        for Te in tes
            Rth = RAH.(CONV * s .* icB[okB], Te)
            cL = Dict(n => sqrt(mean((Rth[lob[okB] .== n] .- r0[okB][lob[okB] .== n]) .^ 2)) for n in lobes)
            c = sum(values(cL)) + λ0L * c0
            if c < best.cost
                best = (cost = c, s, Te, c0, cL, off = Dict(n => dip_offset(Bt[okB], Rth, n) for n in lobes))
            end
        end
    end
    return best
end

scores = [score(Ic[i, :]) for i in eachindex(cands)]
order = sortperm([s.cost for s in scores])
Rav(c) = floor(GEOMETRY[DEV].R) - c.w / 2
TAGN = TARGET === nothing ? "$(NAME)_$(MODE)" : "$(NAME)_$(MODE)_$(TARGET)"
tagbase = joinpath(outdir, TAGN)
lh = join(["rmsR_$(n)L  " for n in lobes])
oh = join(["off$(n) " * " "^4 for n in lobes])
open(tagbase * "_ranking.txt", "w") do io
    println(io, "devices-B/C fit, scan $NAME: $(length(cands)) candidates of device $DEV, target $TGT, physical g, α = 0, kernel T = $TK mK, T_e $MODE (target T_e $(round(TE0; digits = 0)) mK)")
    println(io, "cost = Σ_n≥1 rmsR_nL + $λ0L·rel0L; s = I_c scale, T_N ≈ 0.24 s; offn = lobe-n R-dip centroid − centre (T), + = right")
    println(io, "data: " * join([@sprintf("off%d %+.4f T", n, off_data[n]) for n in lobes], ", "))
    println(io, "rank  g     w   R_av  μ     τΓ    τ      cost    $(lh)rel0L   $(oh)s      T_N≈   T_e")
    for (r, i) in enumerate(order)
        c, s = cands[i], scores[i]
        isfinite(s.cost) || continue
        @printf(io, "%-5d %-5.1f %-3.0f %-5.1f %-5.1f %-5.2f %-6.4f %-7.4f ", r, c.g, c.w, Rav(c), c.μ, c.τΓ, τs[i], s.cost)
        foreach(n -> @printf(io, "%-9.4f ", s.cL[n]), lobes)
        @printf(io, "%-7.4f ", s.c0)
        foreach(n -> @printf(io, "%+-9.4f ", get(s.off, n, NaN)), lobes)
        @printf(io, "%-6.3f %-6.3f %.0f%s\n", s.s, 0.24 * s.s, s.Te, haskey(c, :ξd) ? @sprintf("  ξd %.0f", c.ξd) : "")
    end
end

open(tagbase * "_per_g.txt", "w") do io
    println(io, "best candidate per g (scan $NAME, device $DEV, T_e $MODE); data " * join([@sprintf("off%d %+.4f", n, off_data[n]) for n in lobes], ", "))
    println(io, "g     | w   R_av  μ     τΓ    cost    $(lh)rel0L   $(oh)T_N≈   T_e")
    for gv in sort(unique(getfield.(cands, :g)))
        idx = [i for i in eachindex(cands) if cands[i].g == gv && isfinite(scores[i].cost)]
        isempty(idx) && continue
        i = idx[argmin([scores[j].cost for j in idx])]
        c, s = cands[i], scores[i]
        @printf(io, "%-5.1f | %-3.0f %-5.1f %-5.1f %-5.2f %-7.4f ", gv, c.w, Rav(c), c.μ, c.τΓ, s.cost)
        foreach(n -> @printf(io, "%-9.4f ", s.cL[n]), lobes)
        @printf(io, "%-7.4f ", s.c0)
        foreach(n -> @printf(io, "%+-9.4f ", get(s.off, n, NaN)), lobes)
        @printf(io, "%-6.3f %.0f%s\n", 0.24 * s.s, s.Te, haskey(c, :ξd) ? @sprintf("  ξd %.0f", c.ξd) : "")
    end
end

jldsave(tagbase * "_scores.jld2"; cands, τs, costs = [s.cost for s in scores], s = [s.s for s in scores],
    Te = [s.Te for s in scores], c0 = [s.c0 for s in scores],
    cL = [[get(s.cL, n, NaN) for n in lobes] for s in scores], off = [[get(s.off, n, NaN) for n in lobes] for s in scores],
    lobes, off_data = [off_data[n] for n in lobes], order, dev = DEV, TK, mode = MODE, files, target = TGT)

# contact sheet: best 12
function sheet(idx, file)
    fig = Figure(size = (1500, 280 * cld(length(idx), 3)), fontsize = 12)
    for (k, i) in enumerate(idx)
        c, s = cands[i], scores[i]
        row, col = fldmod1(k, 3)
        g = fig[row, col] = GridLayout()
        ax = Axis(g[1, 1]; ylabel = "I (nA)", titlesize = 11,
            title = @sprintf("#%d g=%.0f w=%.0f (R_av %.1f) μ=%.0f τΓ=%.1f%s T_N≈%.2f T_e=%.0f cost %.3f", k, c.g, c.w, Rav(c), c.μ, c.τΓ,
                haskey(c, :ξd) ? @sprintf(" ξd=%.0f", c.ξd) : "", 0.24 * s.s, s.Te, s.cost))
        scatter!(ax, tg["B"], tg["Isw"]; color = :gray, markersize = 4)
        mm = tg["B"] .> edges[1]
        scatter!(ax, tg["B"][mm], tg["IcAH"][mm] .* (s.Te / tg["Te"]); color = :limegreen, markersize = 5)
        lines!(ax, Bs, CONV * s.s .* Ic[i, :]; color = :orange, linewidth = 2)
        vlines!(ax, edges; color = :gray, linestyle = :dash)
        ylims!(ax, 0, 1.2 * Isw00); xlims!(ax, 0, Bmax)
        axR = Axis(g[2, 1]; xlabel = "B (T)", ylabel = "R0/RN")
        scatter!(axR, Bt, r0; color = :black, markersize = 4)
        lines!(axR, Bs, RAH.(CONV * s.s .* Ic[i, :], s.Te); color = :orchid, linewidth = 2)
        vlines!(axR, edges; color = :gray, linestyle = :dash)
        ylims!(axR, 0, 1.1); xlims!(axR, 0, Bmax)
        rowsize!(g, 2, Relative(0.4)); hidexdecorations!(ax; grid = false, ticks = false)
    end
    save(file, fig)
end
nfin = count(isfinite, [s.cost for s in scores])
if nfin > 0
    sheet(order[1:min(12, nfin)], joinpath("refit_BC", "figures", "$(TAGN)_top.png"))
    # best candidate at each g: how much the fit needs the Zeeman term
    bestg = Int[]
    for gv in sort(unique(getfield.(cands, :g)))
        idx = [i for i in eachindex(cands) if cands[i].g == gv && isfinite(scores[i].cost)]
        isempty(idx) || push!(bestg, idx[argmin([scores[j].cost for j in idx])])
    end
    sheet(bestg, joinpath("refit_BC", "figures", "$(TAGN)_per_g.png"))
else
    println("no candidate covers the target fields (smoke run?): no contact sheet")
end
println(read(tagbase * "_ranking.txt", String)[1:min(end, 3000)])
println(read(tagbase * "_per_g.txt", String))
