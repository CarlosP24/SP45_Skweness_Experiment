function plot_LDOS(ax, name::String; basepath = "data/LDOS")
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack LDOS, params = res
    @unpack Brng, ωrng = params
    ωrng = real.(ωrng)
    ωrng = vcat(ωrng, -reverse(ωrng)[2:end])
    LDOS = cat(LDOS, reverse(LDOS, dims = 2)[:, 2:end], dims = 2)
    heatmap!(ax, Brng, ωrng, LDOS; colormap = :thermal)
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
plot_LDOS(ax, "base_wire")
fig