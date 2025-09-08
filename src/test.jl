function plot_test(da, db; xrng = range(1, 4, 200), φrng = range(0, 2π, 98))
    fig = Figure()

    integranda = mapreduce(permutedims, vcat,  da.(xrng))

    ax = Axis(fig[1, 1]; xlabel = L"x", ylabel = L"\varphi")
    heatmap!(ax, xrng, φrng, abs.(integranda); colormap = :thermal, colorrange = (0, 1e-1))

    ax = Axis(fig[2, 1]; xlabel = L"x", ylabel = L"\varphi")
    heatmap!(ax, xrng, φrng, angle.(integranda); colormap = :redsblues)

    integrandb = mapreduce(permutedims, vcat,  db.(xrng))

    ax = Axis(fig[1, 2]; xlabel = L"x", ylabel = L"\varphi")
    heatmap!(ax, xrng, φrng, abs.(integrandb); colormap = :thermal, colorrange = (0, 1e-1))

    ax = Axis(fig[2, 2]; xlabel = L"x", ylabel = L"\varphi")
    heatmap!(ax, xrng, φrng, angle.(integrandb); colormap = :redsblues)

    return fig
end

plot_test(da, db)