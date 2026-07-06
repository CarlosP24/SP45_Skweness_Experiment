function fig_LDOS_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Combined top two rows from the original figure (height ~188 px, width 600 px)
    fig = Figure(size = (600, 220), fontsize = 16)

    ax = Axis(fig[1, 1];)
    panel_LDOS(ax, name, "left"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, 0.2; text = "Left", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)")
    panel_LDOS(ax, name, "right"; colorrange = (0, 2e-3))
    #hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, -0.2; text = "Right", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[1:2, 1]; ylabel = L"$\\omega$ (meV)", ylabelpadding = -20)
    hidexdecorations!(ax)
    ylims!(ax, (-0.3, 0.3))
    ax.yticks = [-0.2, 0.2]
    hidespines!(ax)

    Colorbar(fig[1:2, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    rowsize!(fig.layout, 1, Relative(1 / 2))
    rowsize!(fig.layout, 2, Relative(1 / 2))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    #rowgap!(fig.layout, 2, 5)

    #style = (font = "CMU Serif Bold", fontsize   = 20)
    #Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -30, 0); style...)

    return fig
end

function fig_resistance_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Middle row from original figure (height ~188 px, width 600 px)
    fig = Figure(size = (600, 188), fontsize = 16)

    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_0/R_N$"
    )
    panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo;)
    xlims!(ax, (first(Trng), last(Trng)))
    hidexdecorations!(ax, ticks = false, grid = false)
    vlines!(ax, [0.11, 0.14, 0.17, 0.195, 0.24]; color = :red, linestyle = :dash)

    rowsize!(fig.layout, 1, Relative(1))
    rowgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "c", padding = (-40, 0, -30, 0); style...)

    return fig
end

function fig_dVdI_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Bottom row from original figure (height ~375 px, width 600 px)
    fig = Figure(size = (600, 375), fontsize = 16)

    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    lims, hmap = panel_dVdI(ax, dVdImat, Trng, Is, Ic, Brng, Ict;)

    Colorbar(fig[1, 2], hmap; label = L"$dV/dI$ (k$\\Omega$)", labelpadding = -20, ticks = ([1000, 1500], ["1", "1.5"]))

    rowsize!(fig.layout, 1, Relative(1))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    Label(fig[1, 1, TopLeft()], "d", padding = (-40, 0, -30, 0); style...)

    return fig
end

function fig_rs(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R/R_N$"
    )
    axd = Axis(fig[2, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$d (R/R_N)/dB$"
    )
    for (i, col) in enumerate(eachcol(rnorm))
        # if mean(col) > 1
        #     continue
        # end
        rs = map(enumerate(col)) do (j, r)
            abs(Is[i]) > IcRs[j] ? NaN : r
        end
        lines!(ax, Trng_Rs, rs; color = (:blue, 0.1))
        abs(Is[i]) < 1e-1 && lines!(ax, Trng_Rs, rs; color = (:red, 0.5), linewidth = 8)

        drs = diff(rs) ./ diff(Trng_Rs)
        lines!(axd, Trng_Rs[2:end], drs .|> abs; color = (:blue, 0.1))

        abs(Is[i]) < 1e-1 && lines!(axd, Trng_Rs[2:end], abs.(drs); color = (:red, 0.5), linewidth = 8)
    end
    hidexdecorations!(ax, ticks = false, grid = false)
    return fig
end
