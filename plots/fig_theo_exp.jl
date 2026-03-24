function ics(name::String; basepath = "data/Js", conv = 74.4)
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec

    return Brng, Ic * conv
end

function process_data(exp_data::Exp_Data, name::String; T::Real = 550, Imax = 120, colorrange = (2.9, 3.2), Blims = (0.08, 0.4), pref = 0.8)
    # Process experimental data 
    data = read_data(exp_data)
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, Imax)
    
    Isw, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Isw)

    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]
    
    i0 = findmin(abs.(Is))[2]

    rnorm = dVdImat ./ RNs
    r0 = rnorm[:, i0] |> vec
    γs = γAH.(r0)
    IcAH = 0.0209 * T .* γs  # T in mK

    Ic = vcat(Isw[1:ilim-1], IcAH[ilim:imax])
    Trng = Ts[1:imax]

    Rs = r0[ilim:imax]
    Trng_Rs = Ts[ilim:imax]

    # Process theory results
    Brng, Ict = ics(name)
    Ict .*= pref
    jmin = findmin(abs.(Brng .- Blims[1]))[2]
    jmax = findmin(abs.(Brng .- Blims[2]))[2]

    R(γ) = besseli0(γ/2)^(-2)
    Rs_theo = R.(Ict[jmin:jmax] ./ (0.0209 * T))

    Brng_Rs = Brng[jmin:jmax]

    return Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm[ilim:imax, :], Ic[ilim:imax], dVdImat[1:imax, :]
end

function panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo; kw...)
    scatter!(ax, Trng_Rs, Rs; color = :black, label = "Experiment", kw...)
    lines!(ax, Brng_Rs, Rs_theo; color = :orange,  label = "Theory", linewidth = 3, kw...)
    axislegend(ax; position = :rb)
end

function panel_dVdI(ax, dVdImat, Ts, Is,  Ic, Brng, Ict; kw... )
    hmap = heatmap!(ax, Ts, Is, abs.(dVdImat); rasterize = 5, colormap = :davos, colorrange = (10^2.9, 10^3.2), colorscale = Makie.pseudolog10, kw...)
    scatter!(ax, Ts, Ic; color = :black, markersize = 8, label = "Experiment")
    scatter!(ax, Ts, -Ic; color = :black, markersize = 8)

    lines!(ax, Brng, Ict; color = :orange, linewidth = 3, label = "Theory")
    lines!(ax, Brng, -Ict; color = :orange, linewidth = 3)
    lims = (first(Ts), last(Ts))
    xlims!(ax, lims)
    axislegend(ax; position = :rt)

    return lims, hmap
end

function panel_LDOS(ax, name::String, side::String; basepath = "data/LDOS", colorrange = (0, 1e-1), kw... )
    # Load resultsº
    path = joinpath(basepath, "fitted_$(name)_$(side).jld2")
    res = load(path)["res"]

    @unpack LDOS, params = res
    @unpack Brng, ωrng = params
    ωrng = real.(ωrng)
    #ωrng = vcat(ωrng, -reverse(ωrng)[2:end])
    #LDOS = cat(LDOSR, reverse(LDOSL, dims = 2)[:, 2:end], dims = 2)
    if side == "left"
        reverse!(LDOS, dims = 2)
        reverse!(ωrng)
        ωrng .*= -1
    end
    heatmap!(ax, Brng, ωrng, LDOS; colormap = :thermal, colorrange, rasterize = 5, kw...)
    hlines!(ax, [0]; color = :white, linewidth = 2)
end


function fig_theo_exp(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    fig = Figure(size = (600, 250 * 3), fontsize = 16,)

    ax = Axis(fig[1, 1]; )
    panel_LDOS(ax, name, "left"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, 0.2; text = "Left", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[2, 1]; )
    panel_LDOS(ax, name, "right"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, -0.2; text = "Right", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[1:2, 1]; ylabel = L"$\omega$ (meV)", ylabelpadding = -20)
    hidexdecorations!(ax)
    ylims!(ax, (-0.3, 0.3))
    ax.yticks = [-0.2,0.2]
    hidespines!(ax)

    Colorbar(fig[1:2, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    ax = Axis(fig[4, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
        )
    lims, hmap =  panel_dVdI(ax, dVdImat, Trng, Is, Ic, Brng, Ict; )

    Colorbar(fig[4, 2], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20, ticks = ([1000, 1500], ["1", "1.5"]))
    ax = Axis(fig[3, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_0/R_N$"
        )
    panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo;)
    xlims!(ax, lims)
    hidexdecorations!(ax, ticks = false, grid = false)
    #vlines!(ax,[0.11, 0.14, 0.17, 0.195, 0.24]; color = :red, linestyle = :dash)

       
    rowsize!(fig.layout, 1, Relative(0.25/2))
    rowsize!(fig.layout, 2, Relative(0.25/2))
    rowsize!(fig.layout, 3, Relative(0.25))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 3, 5)

    style = (font = "CMU Serif Bold", fontsize   = 20)
    Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -30, 0); style...)
    Label(fig[2, 1, TopLeft()], "b",  padding = (-40, 0, -30, 0); style...)
    Label(fig[3, 1, TopLeft()], "c",  padding = (-40, 0, -30, 0); style...)
    Label(fig[4, 1, TopLeft()], "d",  padding = (-40, 0, -30, 0); style...)


    return fig
end

# --- Individual panel figures (sizes match the original `fig_theo_exp` layout) ---
function fig_LDOS_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Combined top two rows from the original figure (height ≃ 188 px, width 600 px)
    fig = Figure(size = (600, 220), fontsize = 16)

    ax = Axis(fig[1, 1]; )
    panel_LDOS(ax, name, "left"; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, 0.2; text = "Left", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[2, 1]; xlabel = L"$B$ (T)")
    panel_LDOS(ax, name, "right"; colorrange = (0, 2e-3))
    #hidexdecorations!(ax, ticks = false, grid = false)
    hideydecorations!(ax)
    text!(ax, 0.36, -0.2; text = "Right", color = :white, align = (:center, :center), fontsize = 20)

    ax = Axis(fig[1:2, 1]; ylabel = L"$\omega$ (meV)", ylabelpadding = -20)
    hidexdecorations!(ax)
    ylims!(ax, (-0.3, 0.3))
    ax.yticks = [-0.2,0.2]
    hidespines!(ax)

    Colorbar(fig[1:2, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    rowsize!(fig.layout, 1, Relative(1/2))
    rowsize!(fig.layout, 2, Relative(1/2))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    #rowgap!(fig.layout, 2, 5)

    #style = (font = "CMU Serif Bold", fontsize   = 20)
    #Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -30, 0); style...)

    return fig
end

function fig_resistance_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Middle row from original figure (height ≃ 188 px, width 600 px)
    fig = Figure(size = (600, 188), fontsize = 16)

    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_0/R_N$"
        )
    panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo;)
    xlims!(ax, (first(Trng), last(Trng)))
    hidexdecorations!(ax, ticks = false, grid = false)
    vlines!(ax,[0.11, 0.14, 0.17, 0.195, 0.24]; color = :red, linestyle = :dash)

    rowsize!(fig.layout, 1, Relative(1))
    rowgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize   = 20)
    Label(fig[1, 1, TopLeft()], "c",  padding = (-40, 0, -30, 0); style...)

    return fig
end

function fig_dVdI_only(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm, IcRs, dVdImat = process_data(exp_data, "fitted_$(name)")

    # Bottom row from original figure (height ≃ 375 px, width 600 px)
    fig = Figure(size = (600, 375), fontsize = 16)

    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
        )
    lims, hmap =  panel_dVdI(ax, dVdImat, Trng, Is, Ic, Brng, Ict; )

    Colorbar(fig[1, 2], hmap; label = L"$dV/dI$ (k$\\Omega$)", labelpadding = -20, ticks = ([1000, 1500], ["1", "1.5"]))

    rowsize!(fig.layout, 1, Relative(1))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)

    style = (font = "CMU Serif Bold", fontsize   = 20)
    Label(fig[1, 1, TopLeft()], "d",  padding = (-40, 0, -30, 0); style...)

    return fig
end

fig = fig_theo_exp(Exp_Data(), "C3")
#save("plots/figures/fig_theo_exp.pdf", fig)
fig

fig_LDOS_only(Exp_Data(), "C3")

##

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
        lines!(axd, Trng_Rs[2:end],  drs .|> abs; color = (:blue, 0.1))

        abs(Is[i]) < 1e-1 && lines!(axd, Trng_Rs[2:end], abs.(drs); color = (:red, 0.5), linewidth = 8)
    end
    hidexdecorations!(ax, ticks = false,  grid = false)
    return fig 
end

fig = fig_rs(Exp_Data(), "C3")
fig