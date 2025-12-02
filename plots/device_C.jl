function fig_deviceC()
    fig = Figure(size = (600, 250 * 4), fontsize = 16,) 

    # Row 1: LDOS left
    ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
    BsL = plot_LDOS(ax, "fitted_C_left", colorrange = (0, 2e-3))
    vlines!(ax, BsL; color = (:blue, .5), linestyle = :dash, label = "LZ left")
    hidexdecorations!(ax, ticks = false)
    text!(ax, 0.37, 0.25; text = "Left", align = (:center, :center), color = :white, fontsize = 20)
    
    # Row 2: LDOS right
    ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)", ylabel = L"$\omega$ (meV)")
    BsR = plot_LDOS(ax, "fitted_C_right", colorrange = (0, 2e-3))
    vlines!(ax, BsR; color = (:blue, .5), linestyle = :dash, label = "LZ right")
    hidexdecorations!(ax, ticks = false)
    text!(ax, 0.37, 0.25; text = "Right", align = (:center, :center), color = :white, fontsize = 20)
    
    # Row 3: Ic
    ax = Axis(fig[3, 1]; xlabel = L"$B$ (T)", ylabel = L"$I_c$ ($\frac{2e}{h}$ meV)", )
    name = "fitted_C"
    plot_Ic(ax, name;  color = :green,  linewidth = 5, label = L"$T=250$mK") 
    ylims!(ax, (0, 2))
    vlines!(ax, BsL; color = (:blue, .5), linestyle = :dash, label = "LZ left")
    #axislegend(ax; position = :rt, framevisible =
    hidexdecorations!(ax, ticks = false)

    csv_path = "Recovered.csv"
    df = CSV.read(csv_path, DataFrame)

    zvals = reverse(Matrix(df[:, 2:end])', dims = 2)

    crit_index = [findmin(zvals[i, :])[2] for i in 1:size(zvals, 1)]

    ax = Axis(fig[4, 1]; xlabel = L"$B$ (T)", ylabel = L"$I_c$ (nA)", backgroundcolor = (:white, 0))
    xrng = range(0, 0.4, length = size(zvals, 1))
    yrng = range(0, 150, length = size(zvals, 2))
    heatmap!(ax, xrng, yrng, zvals, colormap = :redsblues)
    scatter!(ax, xrng, [yrng[i] for i in crit_index]; color = :red,)
    plot_Ic(ax, name;  conv = 74.4, color = :black,  linewidth = 5, label = L"$T=250$mK") 
    ylims!(ax, (0, 150))

    return fig
end

fig = fig_deviceC()
save("plots/figures/deviceC.pdf", fig)
fig


## 
function fig_currents()
    fig = Figure(size = (600, 250 ), fontsize = 16,) 
    name = "C_left_mu=-10_w=20_dw=0_tau=2.5_g=17"
    Bsws = bsws(name)

    ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel =  L"$I_c$ ($\frac{2e}{h}$ meV)", backgroundcolor = (:white, 0) )
    vlines!(ax, Bsws; color = (:black, .5), linestyle = :dash)
    name = "C_mu=25_w=20_dw=0_tau=100_g=17_kBT=0.022_TN=0.19"
    plot_Ic(ax, name;  color = :red,  linewidth = 4,  label = "Little-Parks") 
    name = "fitted_C"
    plot_Ic(ax, name;  color = :green,  linewidth = 4, label = "LP + Carolis") 
    #hidexdecorations!(ax, ticks = false)
    ylims!(ax, 0, 2.5)
    axislegend(ax; position = (1, 1), framevisible = false, orientation = :horizontal, align = (:center, :center))

    style = (font = "CMU Serif Bold", fontsize   = 20)
    Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -35, 0); style...)
    return fig
end

fig = fig_currents()
save("plots/figures/currents_comparison.pdf", fig)
fig

## Export data
ic_total = ics("fitted_C")
ic_LP = ics("C_mu=25_w=20_dw=0_tau=100_g=17_kBT=0.022_TN=0.19")

# Convert vector of vectors to DataFrame with two columns
# Using descriptive column names with units for sharing
df_total = DataFrame("B (T)" => ic_total[1], "Ic (2e/h meV)" => ic_total[2])
df_LP = DataFrame("B (T)" => ic_LP[1], "Ic (2e/h meV)" => ic_LP[2])

CSV.write("ic_total.csv", df_total)
CSV.write("ic_LP.csv", df_LP)