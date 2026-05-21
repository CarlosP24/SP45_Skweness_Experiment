include("plotters/theo_exp_core.jl")
include("plotters/theo_exp_extras.jl")

function fig_theo_exp(exp_data::Exp_Data, name::String; R = 64)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    Φ0 = 2067.833 #T⋅nm^2
    Bns = [1, 2].* Φ0 / (π * R^2)
    Bhalfs = [0.5, 1.5].* Φ0 / (π * R^2)
    BΦs = vcat(Bns, Bhalfs)

    fig = Figure(size = (600, 750), fontsize = 20)

    ax = Axis(fig[1, 1];)
    panel_LDOS(ax, name, "left"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)

    vlines!(ax, Bns; color = :white, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n + 0.4) * Φ0 / (π * R^2), 0.25; text = "$n"*"L", align = (:center, :center), color = :white)
    end

    text!(ax, 0.01, 0.25; text = "Left", color = :white, align = (:left, :center), fontsize = 20)

    gax = Axis(fig[1, 1]; xaxisposition = :top, xlabel = L"\Phi / \Phi_0")
    xlims!(gax, 0, 0.4)
    gax.xticks = (vcat(Bns, Bhalfs), [L"1", L"2", L"0.5", L"1.5"])
    hidespines!(gax)
    hideydecorations!(gax)
    hidexdecorations!(gax, ticks = false, ticklabels = false, label = false)

    ax = Axis(fig[2, 1];)
    panel_LDOS(ax, name, "right"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
 
    vlines!(ax, Bns; color = :white, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    text!(ax, 0.01, -0.25; text = "Right", color = :white, align = (:left, :center), fontsize = 20)


    ax = Axis(fig[1:2, 1]; ylabel = L"$\omega$ (meV)", ylabelpadding = -35)
    hidexdecorations!(ax)
    ylims!(ax, (-0.3, 0.3))
    ax.yticks = [-0.2, 0.2]
    hidespines!(ax)

    Colorbar(fig[1:2, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    ax = Axis(fig[3, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    lims, hmap = panel_dVdI(ax, dVdImat, Trng, Is, Ic, Brng, Ict; colormap = :seismic, lowclip = cgrad(:seismic)[1])
    ax.ylabelpadding = -18

    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n + 0.4 - (n == 2) * 0.1) * Φ0 / (π * R^2), -110; text = "$n"*"L", align = (:center, :center), color = :white)
    end

    hidexdecorations!(ax, ticks = false, grid = false)


    Colorbar(fig[3, 2], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20, ticks = ([1000, 1500], ["1", "1.5"]))
    
    ax = Axis(fig[4, 1];
        xlabel = L"$B$ (T)",        
        ylabel = L"$R_0/R_N$"
    )
    panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo;)
    xlims!(ax, lims)
    
    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n - 0.4) * Φ0 / (π * R^2), 0.5; text = "$n"*"L", align = (:center, :center), color = :black)
    end


    #vlines!(ax,[0.11, 0.14, 0.17, 0.195, 0.24]; color = :red, linestyle = :dash)

    rowsize!(fig.layout, 1, Relative(0.25 / 2))
    rowsize!(fig.layout, 2, Relative(0.25 / 2))
    rowsize!(fig.layout, 4, Relative(0.25 ))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 3, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "a", padding = (-40, 0, -30, 0); style...)
    Label(fig[2, 1, TopLeft()], "b", padding = (-40, 0, -30, 0); style...)
    Label(fig[3, 1, TopLeft()], "c", padding = (-40, 0, -30, 0); style...)
    Label(fig[4, 1, TopLeft()], "d", padding = (-40, 0, -30, 0); style...)

    return fig
end

fig = fig_theo_exp(Exp_Data(), "C3")
save("plots/figures/fig_theo_exp.pdf", fig)
fig
