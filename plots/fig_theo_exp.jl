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
    R0 = dVdImat[:, i0]
    
    Rnorm = R0 ./ RNs
    γs = γAH.(Rnorm)
    IcAH = 0.0209 * T .* γs  # T in mK

    Ic = vcat(Isw[1:ilim-1], IcAH[ilim:imax])
    Trng = Ts[1:imax]

    Rs = Rnorm[ilim:imax]
    Trng_Rs = Ts[ilim:imax]

    # Process theory results
    Brng, Ict = ics(name)
    Ict .*= pref
    jmin = findmin(abs.(Brng .- Blims[1]))[2]
    jmax = findmin(abs.(Brng .- Blims[2]))[2]

    R(γ) = besseli0(γ/2)^(-2)
    Rs_theo = R.(Ict[jmin:jmax] ./ (0.0209 * T))

    Brng_Rs = Brng[jmin:jmax]

    return Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, dVdImat[1:imax, :]
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

function panel_LDOS(ax, name::String; basepath = "data/LDOS", colorrange = (0, 1e-1), kw... )
    # Load resultsº
    pathL = joinpath(basepath, "fitted_$(name)_left.jld2")
    resL = load(pathL)["res"]
    pathR = joinpath(basepath, "fitted_$(name)_right.jld2")
    resR = load(pathR)["res"]

    @unpack LDOS, params = resL
    LDOSL = LDOS
    @unpack LDOS = resR
    LDOSR = LDOS
    @unpack Brng, ωrng = params
    ωrng = real.(ωrng)
    ωrng = vcat(ωrng, -reverse(ωrng)[2:end])
    LDOS = cat(LDOSR, reverse(LDOSL, dims = 2)[:, 2:end], dims = 2)
    heatmap!(ax, Brng, ωrng, LDOS; colormap = :thermal, colorrange, rasterize = 5, kw...)
    hlines!(ax, [0]; color = :white, linewidth = 2)
end


function fig_theo_exp(exp_data::Exp_Data, name::String)
    Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, dVdImat = process_data(Exp_Data(), "fitted_$(name)")

    fig = Figure(size = (600, 250 * 3), fontsize = 16,)

    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$\omega $ (meV)",
        )
    panel_LDOS(ax, name; colorrange = (0, 2e-3))
    hidexdecorations!(ax, ticks = false, grid = false)
    text!(ax, 0.36, 0.2; text = "Left", color = :white, align = (:center, :center), fontsize = 20)
    text!(ax, 0.36, -0.2; text = "Right", color = :white, align = (:center, :center), fontsize = 20)

    Colorbar(fig[1, 2], colormap = :thermal, label = L"LDOS (a. u.) $$", width = 15, limits = (0, 1), ticks = [0, 1], labelpadding = -10)

    ax = Axis(fig[3, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
        )
    lims, hmap =  panel_dVdI(ax, dVdImat, Trng, Is, Ic, Brng, Ict; )

    Colorbar(fig[3, 2], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20, ticks = ([1000, 1500], ["1", "1.5"]))
    ax = Axis(fig[2, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_0/R_N$"
        )
    panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo;)
    xlims!(ax, lims)
    hidexdecorations!(ax, ticks = false, grid = false)

    rowsize!(fig.layout, 1, Relative(0.25))
    rowsize!(fig.layout, 2, Relative(0.25))
    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 2, 5)
    
    style = (font = "CMU Serif Bold", fontsize   = 20)
    Label(fig[1, 1, TopLeft()], "a",  padding = (-40, 0, -30, 0); style...)
    Label(fig[2, 1, TopLeft()], "b",  padding = (-40, 0, -30, 0); style...)
    Label(fig[3, 1, TopLeft()], "c",  padding = (-40, 0, -30, 0); style...)

    return fig
end

fig = fig_theo_exp(Exp_Data(), "C3")
save("plots/figures/fig_theo_exp.pdf", fig)
fig
