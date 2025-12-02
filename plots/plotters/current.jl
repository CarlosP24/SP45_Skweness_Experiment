function plot_Ic(ax, name::String; conv = 1, basepath = "data/Js", kw...)
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    lines!(ax, Brng, conv * Ic; kw...)
    xlims!(ax, (first(Brng), last(Brng)))
    ylims!(ax, (conv * minimum(filter(!isnan, Ic)), conv * maximum(filter(!isnan, Ic))))
end

function ics(name::String; basepath = "data/Js")
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    return [Brng, Ic]
end