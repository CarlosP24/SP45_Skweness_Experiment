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

function plot_LDOS0(ax, name::String; basepath = "data/LDOS", kw... )
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack LDOS, params = res
    @unpack Brng= params

    lines!(ax, Brng, LDOS[:, end]; kw... )
        Φrng = get_Φ(res.wire).(Brng .* cos(res.wire.θ))
    halfints = [h for h in collect(ceil(Int, 2minimum(Φrng))/2 : 0.5 : floor(Int, 2maximum(Φrng))/2) if h != round(h)]
    isws = [argmin(abs.(Φrng .- h)) for h in halfints]
    Bsws = Brng[isws]
    return Bsws
end

function bsws(name::String; basepath = "data/LDOS")
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack params = res
    @unpack Brng = params

    Φrng = get_Φ(res.wire).(Brng .* cos(res.wire.θ))
    halfints = [h for h in collect(ceil(Int, 2minimum(Φrng))/2 : 0.5 : floor(Int, 2maximum(Φrng))/2) if h != round(h)]
    isws = [argmin(abs.(Φrng .- h)) for h in halfints]
    Bsws = Brng[isws]
    return Bsws
end
##
fig = Figure() 
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = "LDOS")
#plot_LDOS(ax, name; colorrange = (0, 2e-3))
plot_LDOS0(ax, "C_left_mu=25_w=20_dw=0_tau=2.5_g=17_T=550"; color = :green, linestyle = :dash)

name = "C_left_mu=20_w=25_tau=2.5_g=17_T=550"
Bsws = plot_LDOS0(ax, name; )
vlines!(ax, Bsws; color = :red, linestyle = :dash)
fig