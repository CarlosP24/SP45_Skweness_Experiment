function fig_conductance(name::String; dir = "data/Conductance")
    data = load(joinpath(dir, name * ".jld2"))["res"]
    @unpack Gs, params = data

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$G$",
        )
    lines!(ax, params.Brng, real.(Gs); color = :blue,)
    lines!(ax, params.Brng, imag.(Gs); color = :red,)

    return fig
end

fig = fig_conductance("fitted_C_test")
fig