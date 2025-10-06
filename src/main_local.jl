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
using CairoMakie
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

wire_systems = Dict(
    [
        key => Wire_System(
            wire = value, 
            params = params_LDOS
        )
        for (key, value) in wires
    ]
)

µ = 10
w = 20
num_mJ = 4

params = Calc_Params(;
    Brng = range(0, 0.4, 101),
    ωrng = range(-0.3, 0, 101) .+ 1e-3im,
)

nameL = "C_ltest_mu=$(μ)_w=$(w)"
nameR = "C_rtest_mu=$(μ)_w=$(w)"
wire_tests = Dict(
    [nameL => Wire_System(
        wire = Params_mm(wires["C_left"];
            μ = µ,
            w = w,
            num_mJ = num_mJ,
            τΓ = 3,
            g = 10,
            α = 0,
        ),
        params = params
    ),
    nameR => Wire_System(
        wire = Params_mm(wires["C_right"];
            μ = µ,
            w = w,
            num_mJ = num_mJ,
            τΓ = 3,
            g = 10,
            α = 0,
        ),
        params = params
    )]
)

merge!(wire_systems, wire_tests)

res = calc_LDOS(nameL)
save(res.path, "res", res)
# res = calc_LDOS(nameR)
# save(res.path, "res", res)

##
fig = Figure()

ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
plot_LDOS(ax, nameL, colorrange = (0, 1e-2))
#hidexdecorations!(ax, ticks = false)

# ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
# plot_LDOS(ax, nameR, colorrange = (0, 2e-2))

fig

##
function plot_Ic(ax, name::String; basepath = "data/Js", kw...)
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    lines!(ax, Brng, Ic; kw...)

end

csv_path = "../Recovered.csv"
df = CSV.read(csv_path, DataFrame)

zvals = reverse( Matrix(df[:, 2:end])', dims = 2)

crit_index = [findmin(zvals[i, :])[2] for i in 1:size(zvals, 1)]

fig = Figure(size = (800, 280))
ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"I_c", backgroundcolor = (:white, 0))
xrng = range(0, 0.4, length = size(zvals, 1))
yrng = range(0, 0.5, length = size(zvals, 2))
heatmap!(ax, xrng, yrng, zvals, colormap = :redsblues)
scatter!(ax, xrng, [yrng[i] for i in crit_index]; color = :red,)

plot_Ic(ax, "C_mu=10_w=15_tau=1.0_g=0"; basepath = "../data/Js", color = :black, linewidth = 5)
xlims!(ax, (0, 0.4))
ylims!(ax, (first(yrng), last(yrng)))
save("Ic_test.pdf", fig)
fig
## 
