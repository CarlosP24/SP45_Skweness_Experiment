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