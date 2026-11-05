# Side-by-side contact sheet of hand-picked candidates from the scans, same scoring as analyze.jl
# (lobe 2 ×2; best overall scale s per candidate).
#   julia --project=plots wip/refit_zeeman/compare.jl  ->  wip/refit_zeeman/compare.png
using JLD2, CairoMakie, Bessels, Statistics, Printf

const CONV, TE = 74.4, 460.0
const L1, L2, λL2, λ0L = (0.085, 0.24), (0.245, 0.38), 2.0, 0.5
const BC2 = 2 * 2067.833848 / (π * 64.0^2)
const S_GRID = exp.(range(log(0.2), log(3.0), 400))
RAH(I) = besseli0(I / (0.0209 * TE) / 2)^(-2)
function linterp(x, y, xq)
    i = searchsortedlast(x, xq); (i < 1 || i >= length(x)) && return NaN
    t = (xq - x[i]) / (x[i+1] - x[i]); return (1 - t) * y[i] + t * y[i+1]
end
function dip_offset(B, R)
    m = [L2[1] <= b <= L2[2] && isfinite(r) for (b, r) in zip(B, R)]
    w = max.(1 .- R[m], 0); sum(w) < 1e-6 && return NaN
    return sum(B[m] .* w) / sum(w) - BC2
end

out = joinpath("wip", "refit_zeeman", "out")
function load_scan(name)
    files = sort([joinpath(out, f) for f in readdir(out) if startswith(f, name * "_") && endswith(f, ".jld2")])
    ps = load.(files)
    return (cands = reduce(vcat, [p["cands"] for p in ps]), Ic = reduce(vcat, [p["Ic"] for p in ps]), Bs = ps[1]["Bs"])
end
scans = Dict(n => load_scan(n) for n in ("scan1", "scan2", "gsweep", "alpha"))
const SET = isempty(ARGS) ? "v1" : ARGS[1]
tg = load(joinpath("wip", "refit_zeeman", "target.jld2"))
Bt = [b for b in tg["B"] if b <= 0.38]; r0 = tg["r0"][[findfirst(==(b), tg["B"]) for b in Bt]]
B0L = [b for (b, I) in zip(tg["Ball"], tg["Iswall"]) if 0.004 <= b < 0.075 && isfinite(I)]
I0L = [I for (b, I) in zip(tg["Ball"], tg["Iswall"]) if 0.004 <= b < 0.075 && isfinite(I)]

picks_alpha = [
    ("α = 0, g = 20 (best lobe 1)", "scan1", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 20.0, α = 0.0)),
    ("α = −10, g = 20 (best overall)", "alpha", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 20.0, α = -10.0)),
    ("α = −30, g = 17.5", "alpha", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 17.5, α = -30.0)),
    ("α = +20, g = 17.5", "alpha", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 17.5, α = 20.0)),
    ("α = 0, g = 48, τΓ = 2.5", "gsweep", (μ = 20.0, w = 25.0, τΓ = 2.5, g = 48.0, α = 0.0)),
    ("α = 0, g = 50, τΓ = 2.5", "gsweep", (μ = 20.0, w = 25.0, τΓ = 2.5, g = 50.0, α = 0.0)),
]
picks_v1 = [
    ("best lobe 1, g = 20", "scan1", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 20.0)),
    ("g = 30, fitted radius", "scan2", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 30.0)),
    ("g = 40, fitted radius", "scan2", (μ = 25.0, w = 25.0, τΓ = 3.0, g = 40.0)),
    ("g = 45, τΓ = 2.5", "scan2", (μ = 25.0, w = 25.0, τΓ = 2.5, g = 45.0)),
    ("published parameters, g = 17π", "scan2", (μ = 20.0, w = 25.0, τΓ = 2.5, g = 17π)),
    ("wide tube, g = 25", "scan2", (μ = 10.0, w = 35.0, τΓ = 2.5, g = 25.0)),
]
picks = SET == "alpha" ? picks_alpha : picks_v1
match(x, c) = all(isapprox(getfield(x, f), getfield(c, f); atol = 1e-6) for f in (:μ, :w, :τΓ, :g)) &&
    isapprox(get(x, :α, 0.0), get(c, :α, 0.0); atol = 1e-6)

fig = Figure(size = (1400, 520 * 2))
for (k, (label, sc, c)) in enumerate(picks)
    S = scans[sc]
    i = findfirst(x -> match(x, c), S.cands)
    i === nothing && (@warn "not found" label; continue)
    ic = S.Ic[i, :]
    icB = [linterp(S.Bs, ic, b) for b in Bt]; ic0 = [linterp(S.Bs, ic, b) for b in B0L]
    ok, ok0 = findall(isfinite, icB), findall(isfinite, ic0)
    best = (c = Inf, s = NaN, l1 = NaN, l2 = NaN)
    for s in S_GRID
        R = RAH.(CONV * s .* icB[ok]); Bo = Bt[ok]
        m1 = [L1[1] <= b <= L1[2] for b in Bo]; m2 = [L2[1] <= b <= L2[2] for b in Bo]
        l1 = sqrt(mean((R[m1] .- r0[ok][m1]) .^ 2)); l2 = sqrt(mean((R[m2] .- r0[ok][m2]) .^ 2))
        c0 = sqrt(mean((CONV * s .* ic0[ok0] .- I0L[ok0]) .^ 2)) / maximum(I0L)
        cc = l1 + λL2 * l2 + λ0L * c0
        cc < best.c && (best = (c = cc, s, l1, l2))
    end
    Rfull = RAH.(CONV * best.s .* ic)
    off = dip_offset(S.Bs, Rfull)
    row, col = fldmod1(k, 3)
    g = fig[row, col] = GridLayout()
    ax = Axis(g[1, 1]; ylabel = "I (nA)", titlesize = 12,
        title = @sprintf("%s\nw=%.0f μ=%.0f τΓ=%.1f T_N≈%.2f | L1 %.3f L2 %.3f off2 %+.4f", label, c.w, c.μ, c.τΓ, 0.24best.s, best.l1, best.l2, off))
    scatter!(ax, tg["Ball"], tg["Iswall"]; color = :gray, markersize = 4)
    scatter!(ax, tg["B"], tg["IcAH"]; color = :green, markersize = 5)
    lines!(ax, S.Bs, CONV * best.s .* ic; color = :orange, linewidth = 2)
    vlines!(ax, [0.0804, 0.2411, BC2]; color = :gray, linestyle = :dash, linewidth = 0.8)
    ylims!(ax, 0, 130); xlims!(ax, 0, 0.4); hidexdecorations!(ax; grid = false, ticks = false)
    axR = Axis(g[2, 1]; xlabel = "B (T)", ylabel = "R0/RN")
    scatter!(axR, Bt, r0; color = :black, markersize = 4)
    lines!(axR, S.Bs, Rfull; color = :orchid, linewidth = 2)
    vlines!(axR, [0.0804, 0.2411, BC2]; color = :gray, linestyle = :dash, linewidth = 0.8)
    ylims!(axR, 0, 1.1); xlims!(axR, 0, 0.4); rowsize!(g, 2, Relative(0.45))
end
Label(fig[0, :], "data lobe-2 dip offset −0.0145 T (dashed: lobe edges and lobe-2 centre 0.321 T)"; fontsize = 14)
save(joinpath("wip", "refit_zeeman", SET == "v1" ? "compare.png" : "compare_$SET.png"), fig)
