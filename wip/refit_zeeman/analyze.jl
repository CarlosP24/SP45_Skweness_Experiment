# Zeeman-corrected refit, step 3: score every scanned candidate against the measurement.
#   julia --project=plots wip/refit_zeeman/analyze.jl <scan name>
# For each candidate the overall scale s of I_c (a transparency rescaling, T_N ≈ 0.24 s, D11) is
# fitted, then the cost is
#   rms[R_th(B) − R_0/R_N(B)] over the target fields (0.085–0.38 T, lobes 1–2)
#   + λ · rms[(s I_c(B) − I_sw(B)) / I_sw(0)] over lobe 0 (B < 0.075 T, the dV/dI peak),
# with R_th = I_0(ħ I_c / 2e k_B T_e)^(−2) as published (theo_exp_core.jl, T_e = 460 mK) and
# I_c in nA through the published conversion 74.4. λ = 0.5.
# Outputs in wip/refit_zeeman/: <name>_ranking.txt, <name>_top.png, <name>_maps.png
using JLD2, CairoMakie, Bessels, Statistics, Printf

const NAME = isempty(ARGS) ? "scan1" : ARGS[1]
const CONV = 74.4          # code units -> nA (plots/plotters/theo_exp_core.jl)
const TE = 460.0           # mK
const λ0L = 0.5
const S_GRID = exp.(range(log(0.2), log(3.0), 400))

outdir = joinpath("wip", "refit_zeeman", "out")
files = isfile(joinpath(outdir, "$NAME.jld2")) ? [joinpath(outdir, "$NAME.jld2")] :
    sort([joinpath(outdir, f) for f in readdir(outdir) if startswith(f, NAME * "_") && endswith(f, ".jld2")])
isempty(files) && error("no output for $NAME in $outdir")
parts = load.(files)
d = first(parts)
all(p -> p["Bs"] == d["Bs"], parts) || error("parts of $NAME use different field grids")
cands = reduce(vcat, [p["cands"] for p in parts])
Ic = reduce(vcat, [p["Ic"] for p in parts])
Bs = d["Bs"]
println("merged $(length(files)) file(s), $(length(cands)) candidates")
tg = load(joinpath("wip", "refit_zeeman", "target.jld2"))

"linear interpolation of y(x) at xq, NaN outside or across NaN"
function linterp(x, y, xq)
    i = searchsortedlast(x, xq)
    (i < 1 || i >= length(x)) && return NaN
    t = (xq - x[i]) / (x[i+1] - x[i])
    return (1 - t) * y[i] + t * y[i+1]
end

RAH(Ic_nA) = besseli0(Ic_nA / (0.0209 * TE) / 2)^(-2)

# target
Bt = [b for b in tg["B"] if b <= maximum(Bs)]
it = [findfirst(==(b), tg["B"]) for b in Bt]
r0 = tg["r0"][it]
B0L = [b for (b, I) in zip(tg["Ball"], tg["Iswall"]) if 0.004 <= b < 0.075 && isfinite(I)]
Isw0L = [I for (b, I) in zip(tg["Ball"], tg["Iswall"]) if 0.004 <= b < 0.075 && isfinite(I)]
Isw00 = maximum(Isw0L)

const L1 = (0.085, 0.24)    # lobe 1 in the target window (T)
const L2 = (0.245, 0.38)    # lobe 2
const λL2 = 2.0             # lobe 2 weighted ×2 (Carlos, 2026-11-05)
const BC2 = 2 * 2067.833848 / (π * 64.0^2)   # centre of lobe 2, Φ = 2Φ0 at R_LP = 64 nm (T)

"centroid of the lobe-2 resistance dip, weight max(1 − R, 0), minus the lobe centre (T); < 0 = left skew"
function dip_offset(B, R)
    m = [L2[1] <= b <= L2[2] && isfinite(r) for (b, r) in zip(B, R)]
    w = max.(1 .- R[m], 0)
    sum(w) < 1e-6 && return NaN
    return sum(B[m] .* w) / sum(w) - BC2
end

function score(icurve)
    best = (cost = Inf, s = NaN, cR = NaN, c0 = NaN, cL1 = NaN, cL2 = NaN, off2 = NaN)
    icB = [linterp(Bs, icurve, b) for b in Bt]
    ic0 = [linterp(Bs, icurve, b) for b in B0L]
    okB, ok0 = findall(isfinite, icB), findall(isfinite, ic0)
    # a few failed fields are tolerated; a candidate missing >10% of either set is rejected
    (length(okB) < 0.9 * length(icB) || length(ok0) < 0.9 * length(ic0)) && return best
    for s in S_GRID
        Rth = RAH.(CONV * s .* icB[okB])
        cR = sqrt(mean((Rth .- r0[okB]) .^ 2))
        c0 = sqrt(mean(((CONV * s .* ic0[ok0]) .- Isw0L[ok0]) .^ 2)) / Isw00
        Bok = Bt[okB]
        m1 = [L1[1] <= b <= L1[2] for b in Bok]
        m2 = [L2[1] <= b <= L2[2] for b in Bok]
        cL1 = sqrt(mean((Rth[m1] .- r0[okB][m1]) .^ 2))
        cL2 = sqrt(mean((Rth[m2] .- r0[okB][m2]) .^ 2))
        c = cL1 + λL2 * cL2 + λ0L * c0
        if c < best.cost
            best = (cost = c, s, cR, c0, cL1, cL2, off2 = dip_offset(Bok, Rth))
        end
    end
    return best
end

scores = [score(Ic[i, :]) for i in eachindex(cands)]
order = sortperm([s.cost for s in scores])

off2_data = dip_offset(Bt, r0)
open(joinpath("wip", "refit_zeeman", "$(NAME)_ranking.txt"), "w") do io
    println(io, "Zeeman-corrected refit, scan $NAME: $(length(cands)) candidates, physical g, T = $(d["T"]) mK")
    println(io, "cost = rmsR_L1 + $λL2·rmsR_L2 + $λ0L·rel0L; s = I_c scale, T_N ≈ 0.24 s")
    @printf(io, "lobe-2 dip offset (centroid − lobe centre %.3f T): data %+.4f T (negative = left skew)\n", BC2, off2_data)
    println(io, "rank   g     w(nm) R_av   μ(meV) τΓ    α      cost     rmsR_L1  rmsR_L2  off2     rel0L    s      T_N≈")
    for (r, i) in enumerate(order)
        c, s = cands[i], scores[i]
        @printf(io, "%-6d %-5.1f %-5.0f %-6.1f %-6.1f %-5.2f %-6.1f %-8.4f %-8.4f %-8.4f %-+8.4f %-8.4f %-6.3f %-6.3f\n",
            r, c.g, c.w, 60 - c.w / 2, c.μ, c.τΓ, get(c, :α, 0.0), s.cost, s.cL1, s.cL2, s.off2, s.c0, s.s, 0.24 * s.s)
    end
end

# best per (α, g) when α varies: overall cost, and the best lobe-2 fit keeping lobe 1 (rmsR_L1 ≤ 0.12)
if length(unique(get.(cands, :α, 0.0))) > 1
    open(joinpath("wip", "refit_zeeman", "$(NAME)_per_alpha.txt"), "w") do io
        @printf(io, "best candidate per (α, g) (scan %s); data lobe-2 dip offset %+.4f T\n", NAME, off2_data)
        println(io, "α      g     | best cost: w  μ   τΓ   cost   L1     L2     off2    | best L2 with L1≤0.12: w  μ   τΓ   L1     L2     off2")
        for αv in sort(unique(get.(cands, :α, 0.0))), gv in sort(unique(getfield.(cands, :g)))
            idx = [i for i in eachindex(cands) if get(cands[i], :α, 0.0) == αv && cands[i].g == gv && isfinite(scores[i].cost)]
            isempty(idx) && continue
            a = idx[argmin([scores[i].cost for i in idx])]
            good = [i for i in idx if scores[i].cL1 <= 0.12]
            ca, sa = cands[a], scores[a]
            @printf(io, "%-6.1f %-5.1f | %-3.0f %-4.0f %-4.2f %.4f %.4f %.4f %+.4f | ", αv, gv, ca.w, ca.μ, ca.τΓ, sa.cost, sa.cL1, sa.cL2, sa.off2)
            if isempty(good)
                println(io, "none")
            else
                b = good[argmin([scores[i].cL2 for i in good])]
                cb, sb = cands[b], scores[b]
                @printf(io, "%-3.0f %-4.0f %-4.2f %.4f %.4f %+.4f\n", cb.w, cb.μ, cb.τΓ, sb.cL1, sb.cL2, sb.off2)
            end
        end
    end
end

# best per g: overall cost, and the best lobe-2 fit that keeps lobe 1 (rmsR_L1 ≤ 0.12)
open(joinpath("wip", "refit_zeeman", "$(NAME)_per_g.txt"), "w") do io
    @printf(io, "best candidate per g (scan %s); data lobe-2 dip offset %+.4f T\n", NAME, off2_data)
    println(io, "g      | best cost: w  μ   τΓ   cost   L1     L2     off2    | best L2 with L1≤0.12: w  μ   τΓ   L1     L2     off2")
    for gv in sort(unique(getfield.(cands, :g)))
        idx = [i for i in eachindex(cands) if cands[i].g == gv && isfinite(scores[i].cost)]
        isempty(idx) && continue
        a = idx[argmin([scores[i].cost for i in idx])]
        good = [i for i in idx if scores[i].cL1 <= 0.12]
        ca, sa = cands[a], scores[a]
        @printf(io, "%-6.1f | %-3.0f %-4.0f %-4.2f %.4f %.4f %.4f %+.4f | ", gv, ca.w, ca.μ, ca.τΓ, sa.cost, sa.cL1, sa.cL2, sa.off2)
        if isempty(good)
            println(io, "none")
        else
            b = good[argmin([scores[i].cL2 for i in good])]
            cb, sb = cands[b], scores[b]
            @printf(io, "%-3.0f %-4.0f %-4.2f %.4f %.4f %+.4f\n", cb.w, cb.μ, cb.τΓ, sb.cL1, sb.cL2, sb.off2)
        end
    end
end

# lobe 2 alone (with each candidate's best overall scale)
order2 = sortperm([isfinite(s.cL2) ? s.cL2 : Inf for s in scores])
open(joinpath("wip", "refit_zeeman", "$(NAME)_ranking_L2.txt"), "w") do io
    println(io, "Ranked by rms[R_th − R0/RN] in lobe 2 only ($(L2[1])–$(L2[2]) T); flat R = 1 scores $(round(sqrt(mean((1 .- r0[[L2[1] <= b <= L2[2] for b in Bt]]) .^ 2)); digits = 4))")
    println(io, "rank   g     w(nm) μ(meV) τΓ    rmsR_L2  rmsR_L1  cost     T_N≈")
    for (r, i) in enumerate(order2[1:min(40, end)])
        c, s = cands[i], scores[i]
        @printf(io, "%-6d %-5.1f %-5.0f %-6.1f %-5.2f %-8.4f %-8.4f %-8.4f %-6.3f\n",
            r, c.g, c.w, c.μ, c.τΓ, s.cL2, s.cL1, s.cost, 0.24 * s.s)
    end
end

# contact sheets: best 12 overall, best 9 in lobe 2 alone
function sheet(idx, file)
ntop = length(idx)
fig = Figure(size = (1400, 260 * cld(ntop, 3)))
for (k, i) in enumerate(idx)
    c, s = cands[i], scores[i]
    row, col = fldmod1(k, 3)
    g = fig[row, col] = GridLayout()
    ax = Axis(g[1, 1]; ylabel = "I (nA)", title = @sprintf("#%d  g=%.1f α=%.0f w=%.0f μ=%.0f τΓ=%.2f  T_N≈%.2f  L2 %.3f off2 %+.4f",
        k, c.g, get(c, :α, 0.0), c.w, c.μ, c.τΓ, 0.24 * s.s, s.cL2, s.off2), titlesize = 11)
    scatter!(ax, tg["Ball"], tg["Iswall"]; color = :gray, markersize = 4)
    scatter!(ax, tg["B"], tg["IcAH"]; color = :green, markersize = 5)
    lines!(ax, Bs, CONV * s.s .* Ic[i, :]; color = :orange, linewidth = 2)
    ylims!(ax, 0, 130); xlims!(ax, 0, 0.4)
    axR = Axis(g[2, 1]; xlabel = "B (T)", ylabel = "R0/RN")
    scatter!(axR, Bt, r0; color = :black, markersize = 4)
    lines!(axR, Bs, RAH.(CONV * s.s .* Ic[i, :]); color = :orchid, linewidth = 2)
    ylims!(axR, 0, 1.1); xlims!(axR, 0, 0.4)
    rowsize!(g, 2, Relative(0.4)); hidexdecorations!(ax; grid = false, ticks = false)
end
save(joinpath("wip", "refit_zeeman", file), fig)
end
sheet(order[1:min(12, end)], "$(NAME)_top.png")
sheet(order2[1:min(9, end)], "$(NAME)_top_L2.png")

# cost maps (w, μ) for each g, minimised over τΓ
gs, ws, μs = sort(unique(getfield.(cands, :g))), sort(unique(getfield.(cands, :w))), sort(unique(getfield.(cands, :μ)))
fig2 = Figure(size = (420 * length(gs), 380))
for (k, gv) in enumerate(gs)
    M = fill(NaN, length(ws), length(μs))
    for (i, c) in enumerate(cands)
        c.g == gv || continue
        a, b = findfirst(==(c.w), ws), findfirst(==(c.μ), μs)
        M[a, b] = isnan(M[a, b]) ? scores[i].cost : min(M[a, b], scores[i].cost)
    end
    ax = Axis(fig2[1, k]; xlabel = "w (nm)", ylabel = "μ (meV)", title = "g = $gv (min over τΓ)",
        xticks = ws, yticks = μs)
    fin = filter(isfinite, [s.cost for s in scores])
    hm = heatmap!(ax, ws, μs, M; colormap = :viridis, colorrange = isempty(fin) ? (0, 1) : (minimum(fin), max(minimum(fin) + 0.05, 0.5)))
    k == length(gs) && Colorbar(fig2[1, k+1], hm; label = "cost")
end
save(joinpath("wip", "refit_zeeman", "$(NAME)_maps.png"), fig2)
println(read(joinpath("wip", "refit_zeeman", "$(NAME)_ranking.txt"), String)[1:min(end, 2500)])
