# Report figure: lobe skew of the best candidate at each g against the measured skew, devices B and C.
#   julia --project=plots refit_BC/fig_skew_vs_g.jl [T_e mode]
# Skew = centroid of the lobe's R_0/R_N dip (weight max(1 − R, 0)) minus the lobe centre (analyze.jl).
# Output: refit_BC/figures/skew_vs_g.png
using JLD2, CairoMakie, Printf

const MODE = isempty(ARGS) ? "fixed" : ARGS[1]
fig = Figure(size = (1100, 460), fontsize = 15)
for (k, dev) in enumerate(("B", "C"))
    f = joinpath("refit_BC", "out", "$(dev)_$(MODE)_scores.jld2")
    isfile(f) || continue
    sc = load(f)
    cands, costs, off, lobes = sc["cands"], sc["costs"], sc["off"], sc["lobes"]
    gs = sort(unique(getfield.(cands, :g)))
    ax = Axis(fig[1, k]; xlabel = "physical g", ylabel = "R-dip centroid − lobe centre (mT)",
        title = "device $dev (best candidate at each g, T_e $MODE)", xticks = gs)
    colors = Makie.wong_colors()
    for (j, n) in enumerate(lobes)
        th = Float64[]
        for gv in gs
            idx = [i for i in eachindex(cands) if cands[i].g == gv && isfinite(costs[i])]
            push!(th, isempty(idx) ? NaN : 1e3 * off[idx[argmin(costs[idx])]][j])
        end
        scatterlines!(ax, gs, th; color = colors[j], label = "$(n)L theory", markersize = 10)
        hlines!(ax, [1e3 * sc["off_data"][j]]; color = colors[j], linestyle = :dash, label = "$(n)L data")
    end
    hlines!(ax, [0]; color = :gray, linewidth = 0.5)
    Legend(fig[2, k], ax; orientation = :horizontal, nbanks = 2, labelsize = 12, framevisible = false)
    ax2 = Axis(fig[1, k]; yaxisposition = :right, ylabel = "best cost", ylabelcolor = :gray50, yticklabelcolor = :gray50)
    hidexdecorations!(ax2); hidespines!(ax2)
    best = [minimum(costs[[i for i in eachindex(cands) if cands[i].g == gv]]) for gv in gs]
    lines!(ax2, gs, best; color = :gray50, linestyle = :dot, linewidth = 2)
end
save(joinpath("refit_BC", "figures", "skew_vs_g.png"), fig)
println("saved refit_BC/figures/skew_vs_g.png")
