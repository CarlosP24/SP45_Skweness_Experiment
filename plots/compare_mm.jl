# Check of the m_J-independent builder against the published multimode results.
#   julia --project=plots plots/compare_mm.jl [dir]     (dir with Js/ and LDOS/, default "data")
# Compares <dir>/Js/{fitted,depleted}_C3.jld2 and <dir>/LDOS/fitted_C3_{left,right}.jld2 with the
# multimode results kept as data/Js/*_mm.jld2 and data/LDOS/*_mm.jld2 (build_cyl_mm, "g = 17" in its
# convention, i.e. a physical g = 17π). Prints the worst I_c deviation relative to the peak and LDOS
# norms, and saves figures/compare_mm.png.
using JLD2, CairoMakie, Statistics, Printf

const DIR = isempty(ARGS) ? "data" : ARGS[1]
Ic(res) = vec(maximum(res.Js; dims = 2))
norm2(x) = sqrt(sum(abs2, x))

fig = Figure(size = (1400, 800), fontsize = 14)
for (k, name) in enumerate(("fitted_C3", "depleted_C3"))
    new, old = load(joinpath(DIR, "Js", "$name.jld2"))["res"], load(joinpath("data", "Js", "$(name)_mm.jld2"))["res"]
    Bn, Bo = collect(new.params.Brng), collect(old.params.Brng)
    Bn == Bo || error("$name: field grids differ")
    In, Io = Ic(new), Ic(old)
    ok = isfinite.(In) .& isfinite.(Io)
    dev = abs.(In[ok] .- Io[ok]) ./ maximum(Io[ok])
    @printf("I_c %-12s worst |new − mm| / peak = %.2e (at B = %.4f T), rms %.2e, peak new/mm = %.5f / %.5f\n",
        name, maximum(dev), Bn[ok][argmax(dev)], sqrt(mean(dev .^ 2)), maximum(In[ok]), maximum(Io[ok]))
    ax = Axis(fig[1, k]; title = "$name: I_c", xlabel = "B (T)")
    lines!(ax, Bo, Io; color = :black, linewidth = 3, label = "multimode (published)")
    lines!(ax, Bn, In; color = :orange, linestyle = :dash, linewidth = 2, label = "m_J sum")
    axislegend(ax)
end
for (k, side) in enumerate(("left", "right"))
    new, old = load(joinpath(DIR, "LDOS", "fitted_C3_$side.jld2"))["res"], load(joinpath("data", "LDOS", "fitted_C3_$(side)_mm.jld2"))["res"]
    Ln, Lo = new.LDOS, old.LDOS
    size(Ln) == size(Lo) || error("LDOS $side: sizes differ")
    ok = isfinite.(Ln) .& isfinite.(Lo)      # both are singular at the same gapless points (right lead, B = 0.4 T)
    @printf("LDOS %-5s ‖new − mm‖/‖mm‖ = %.2e, max |new − mm| / max(mm) = %.2e (%d non-finite points, new %d / mm %d)\n", side,
        norm2(Ln[ok] .- Lo[ok]) / norm2(Lo[ok]), maximum(abs.(Ln[ok] .- Lo[ok])) / maximum(Lo[ok]), count(!, ok), count(!isfinite, Ln), count(!isfinite, Lo))
    Lo = ifelse.(ok, Lo, NaN); Ln = ifelse.(ok, Ln, NaN)
    ax = Axis(fig[2, k]; title = "LDOS $side: (new − mm) / max(mm)", xlabel = "B (T)", ylabel = "ω (meV)")
    hm = heatmap!(ax, collect(new.params.Brng), real.(collect(new.params.ωrng)), (Ln .- Lo) ./ maximum(Lo[ok]);
        colormap = :balance, colorrange = (-0.05, 0.05))
    k == 2 && Colorbar(fig[2, 3], hm)
end
mkpath("figures")
save(joinpath("figures", "compare_mm_$(basename(DIR)).png"), fig)
println("saved figures/compare_mm_$(basename(DIR)).png")
