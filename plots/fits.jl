function plot_LDOS(ax, name::String; basepath = "data/LDOS", colorrange = (0, 1e-1))
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack LDOS, params = res
    @unpack Brng, ωrng = params
    ωrng = real.(ωrng)
    ωrng = vcat(ωrng, -reverse(ωrng)[2:end])
    LDOS = cat(LDOS, reverse(LDOS, dims = 2)[:, 2:end], dims = 2)
    heatmap!(ax, Brng, ωrng, LDOS; colormap = :thermal, colorrange)

    Φrng = get_Φ(res.wire).(Brng .* cos(res.wire.θ))
    halfints = [h for h in collect(ceil(Int, 2minimum(Φrng))/2 : 0.5 : floor(Int, 2maximum(Φrng))/2) if h != round(h)]
    isws = [argmin(abs.(Φrng .- h)) for h in halfints]
    Bsws = Brng[isws]
    return Bsws
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
plot_LDOS(ax, "C_right")
fig

##
function plot_Ic(ax, name::String; basepath = "data/Js")
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    lines!(ax, Brng, Ic;)

end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"I_c")
plot_Ic(ax, "C_1.0")
fig

##
fig = Figure(size = (500, 600), fontsize = 20)
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
Bsleft = plot_LDOS(ax, "C_left", colorrange = (0, 4e-4))
vlines!(ax, Bsleft, color = :orange, linestyle = :dash, linewidth = 3)

hidexdecorations!(ax, ticks = false)
ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
Bsright = plot_LDOS(ax, "C_right", colorrange = (0, 5e-5))
vlines!(ax, Bsright, color = :green, linestyle = :dash, linewidth = 3)

hidexdecorations!(ax, ticks = false)
ax = Axis(fig[3, 1]; xlabel = L"$B$ (T)", ylabel = L"I_c")
plot_Ic(ax, "C_1.0")
vlines!(ax, Bsleft; color = :orange, linestyle = :dash, linewidth = 3)
vlines!(ax, Bsright; color = :green, linestyle = :dash, linewidth = 3)
xlims!(ax, (0, 0.6))
fig