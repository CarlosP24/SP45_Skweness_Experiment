using Distributed
addprocs(8)

@everywhere begin
    using Pkg; Pkg.activate(".");
    Pkg.instantiate(); Pkg.precompile()
end 
    
cd("src")

using Pkg, JLD2
using CSV, DataFrames
include("utilities.jl")
ensure_pkgs(["Quantica", "FullShell", "Parameters", "ProgressMeter", "Interpolations", "Roots"])
@everywhere begin
    using Quantica, FullShell, JLD2
    using ProgressMeter, Parameters
    using Interpolations, Roots
    # Load code
    include("utilities.jl")
    include("models/params.jl")
    include("models/wires.jl")
    include("models/junctions.jl")
    include_all("builders")
    include_all("operators")
    include_all("parallelizers")
    include_all("calculations")
end

##
using CairoMakie, JLD2, Parameters
using CSV, DataFrames
cd("src")
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
##
common = "mu=25.0_w=20_dw=30_tau=2.5_g=17"
nameL = "C_left_$(common)"
nameR = "C_right_$(common)"

read_wire(nameL)
read_wire(nameR)

res = calc_LDOS(nameL)
save(res.path, "res", res)
res = calc_LDOS(nameR)
save(res.path, "res", res)

##
fig = Figure(size = (600, 600))

ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
plot_LDOS(ax, nameL, colorrange = (0, 4e-3))
text!(ax, 0.37, 0.25; text = "Left", align = (:center, :center), color = :white, fontsize = 20)
# hidexdecorations!(ax, ticks = false)

ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
plot_LDOS(ax, nameR, colorrange = (0, 4e-3))
text!(ax, 0.37, 0.25; text = "Right", align = (:center, :center), color = :white, fontsize = 20)
save("../plots/study_figs/LDOS_$(common).pdf", fig)
fig

##
function plot_Ic(ax, name::String; conv = 1, basepath = "data/Js", kw...)
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    scatter!(ax, Brng, conv * Ic; kw...)
    xlims!(ax, (first(Brng), last(Brng)))
    ylims!(ax, (0, conv * maximum(filter(!isnan, Ic))))
end

csv_path = "../Recovered.csv"
df = CSV.read(csv_path, DataFrame)

zvals = reverse(Matrix(df[:, 2:end])', dims = 2)

crit_index = [findmin(zvals[i, :])[2] for i in 1:size(zvals, 1)]

fig = Figure(size = (800, 280))
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$I_c$ (nm)", backgroundcolor = (:white, 0))
xrng = range(0, 0.4, length = size(zvals, 1))
yrng = range(0, 150, length = size(zvals, 2))
heatmap!(ax, xrng, yrng, zvals, colormap = :redsblues)
scatter!(ax, xrng, [yrng[i] for i in crit_index]; color = :red,)

name = "C_mu=25_w=20_dw=0_tau=2.5_g=17_kBT=0.022"
plot_Ic(ax, name; conv = 74.4, basepath = "../data/Js", color = :lightgreen,)
#xlims!(ax, (0, 0.4))
#ylims!(ax, (first(yrng), last(yrng)))
ylims!(ax, 0, 150)

save("../plots/study_figs/$(name).pdf", fig)
fig
## 

fig = Figure(size = (800, 280))
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$I_c$ (nA)")
name = "C_mu=25.0_w=20_dw=0_tau=2.5_g=17"
plot_Ic(ax, name; conv = 74.4, basepath = "../data/Js", color = :black, linewidth = 5, label = L"T=0") #  meV * 2e / h = 74.4 nA
name = "C_mu=25_w=20_dw=0_tau=2.5_g=17_kBT=0.022"
plot_Ic(ax, name; conv = 74.4, basepath = "../data/Js", color = :green, linestyle = :dash, linewidth = 5, label = L"$T=250$mK") #  meV * 2e / h = 74.4 nA
ylims!(ax, 0, 75)
axislegend(ax; position = :rt, framevisible = false, textsize = 20)
fig


