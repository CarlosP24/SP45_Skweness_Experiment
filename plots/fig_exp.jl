function fig_exp(exp_dataB::Exp_Data, exp_dataT::Exp_Data; R = 64, Imax = 800, Bmax = 0.4)

    Φ0 = 2067.833 #T⋅nm^2
    Bns = [1, 2].* Φ0 / (π * R^2)
    Bhalfs = [0.5, 1.5, 2.5].* Φ0 / (π * R^2)
    Bs, Is, dVdImat, _, _ = prepare_data(exp_dataB, Imax)
    iB = findmin(abs.(Bs .- Bmax))[2]

    fig = Figure(size = (600, 250 * 2 + 50), fontsize = 20)
    ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$I$ (nA)")
    ax.yticks = -600:300:600
    dVdI = abs.(dVdImat[1:iB, :])

    hmap = heatmap!(ax, Bs[1:iB], Is, dVdI; rasterize = 5, colormap = :seismic, colorrange = (0.4e3, 1.6e3), lowclip = cgrad(:seismic)[1])
    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n + 0.4) * Φ0 / (π * R^2), -650; text = "$n"*"L", align = (:center, :center), color = :white)
    end

    text!(ax, 0.25, 400; text = "Skewed\n switching current", align = (:center, :center), color = :black, fontsize = 16)
    arrows2d!(ax, [0.25], [220], [0.23], [50]; argmode = :endpoint, shaftwidth = 2, tipwidth = 10)
    arrows2d!(ax, [0.25], [220], [0.27], [50]; argmode = :endpoint, shaftwidth = 2, tipwidth = 10)

    text!(ax, 0.1, 460; text = "MAR", align = (:center, :center), color = :black, fontsize = 16)
    arrows2d!(ax, [0.08], [450], [0.055], [310]; argmode = :endpoint, shaftwidth = 2, tipwidth = 10)
    #arrows2d!(ax, [0.1], [450], [0.12], [310]; argmode = :endpoint)


    gax = Axis(fig[1, 1]; xaxisposition = :top, xlabel = L"\Phi / \Phi_0")
    xlims!(gax, Bs[1], Bs[iB])
    gax.xticks = (vcat(Bns, Bhalfs), [L"1", L"2", L"0.5", L"1.5", L"2.5"])
    hidespines!(gax)
    hideydecorations!(gax)
    hidexdecorations!(gax, ticks = false, ticklabels = false, label = false)
    
    Colorbar(fig[1, 2], hmap; ticks = ([500, 1500], ["0.5", "1.5"]), label = L"$dV/dI$ (k$\Omega$)", labelpadding = -15)

    hidexdecorations!(ax, ticks = false)

    Ts, Bs, dVdImat, _, _ = prepare_data(exp_dataT, 0.6; Icoef = 1)
    iB = findmin(abs.(Bs .- Bmax))[2]
    ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)", ylabel = L"$T$ (K)")
    gax = Axis(fig[2, 1]; xaxisposition = :top)
    xlims!(gax, Bs[1], Bs[iB])
    gax.xticks = (vcat(Bns, Bhalfs), [L"1", L"2", L"0.5", L"1.5", L"2.5"])
    hidespines!(gax)
    hideydecorations!(gax)
    hidexdecorations!(gax, ticks = false,)

    dVdI = abs.(dVdImat[:, 1:iB])'

    hmap = heatmap!(ax, Bs[1:iB], Ts, dVdI; rasterize = 5, colormap = :seismic, colorrange = (0.6e3, 1.6e3), lowclip = cgrad(:seismic)[1])
    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n + 0.4) * Φ0 / (π * R^2), 1.25; text = "$n"*"L", align = (:center, :center), color = :white)
    end

    text!(ax, 0.28, 0.45; text = "Skewed\n critical\n current", align = (:center, :center), color = :white, fontsize = 16)
    arrows2d!(ax, [0.255], [0.44], [0.22], [0.4]; argmode = :endpoint, color = :white, shaftwidth = 2, tipwidth = 10)
   
    #text!(ax, 0.088, 1.07; text = L"T_c", align = (:center, :center), color = :white, fontsize = 20)
    #arrows2d!(ax, [0.085], [1.0], [0.1], [0.9]; argmode = :endpoint, color = :white, shaftwidth = 2, tipwidth = 10)
    #arrows2d!(ax, [0.085], [1.0], [0.065], [0.9]; argmode = :endpoint, color = :white, shaftwidth = 2, tipwidth = 10)

    Colorbar(fig[2, 2], hmap; ticks = ([700, 1500], ["0.7", "1.5"]), label = L"$dV/dI$ (k$\Omega$)", labelpadding = -15, )
    rowgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -30, 0); style...)
    Label(fig[2, 1, TopLeft()], "b", padding = (-40, 0, -30, 0); style...)
    return fig
end

fig = fig_exp(Exp_Data(;
    keysight = "data_exp/250mK_DEV13_277_18_dmm_keysight_signal.dat",
    lockin = "data_exp/250mK_DEV13_277_18_mf1_X.dat"
    ),
    Exp_Data(;
    keysight = "data_exp/250mK_DEV13_294_35_V_drop.dat",
    lockin = "data_exp/250mK_DEV13_294_35_R_demod.dat"
    )
)
save("plots/figures/fig_exp.pdf", fig)
fig