# ============================================================================
# Common utility functions
# ============================================================================

function γAH(R)
    """Compute the Ambegaokar-Helkin normalized current γ from normal resistance ratio R."""
    if R < 1
        f(γ) = R - besseli0(γ/2)^(-2)
        γ = find_zero(f, (0, 100);)
        return γ
    end
    return NaN
end

function select_top_peaks(all_indices::Vector, Is::Vector, n::Int=2)
    """Select the n highest peaks by Is value."""
    if length(all_indices) > n
        sorted_idx = sortperm(Is[all_indices], rev=true)
        return all_indices[sorted_idx[1:n]]
    else
        return all_indices
    end
end

function get_peak_indices_smoothed(dVdI_row::Vector; window=101, order=1, min_prom=6)
    """Find peak indices in a smoothed dVdI curve."""
    rf = savitzky_golay(dVdI_row, window, order).y
    pks = findmaxima(rf;)
    pks = peakproms(pks; min = min_prom)
    return pks.indices
end

function extract_peak_currents(dVdImat::Matrix, Is::Vector, Ts::Vector)
    """Extract critical currents from peak analysis across all B field values."""
    Ics = map(axes(dVdImat, 1)) do ir
        all_indices = get_peak_indices_smoothed(dVdImat[ir, :])
        selected = select_top_peaks(all_indices, Is, 2)
        
        if length(selected) == 2
            Ic_u, Ic_l = Is[selected[1]], Is[selected[2]]
            if Ic_u < Ic_l
                Ic_u, Ic_l = Ic_l, Ic_u
            end
            return Ic_u, Ic_l
        elseif length(selected) == 1
            if Is[selected[1]] > 0
                return Is[selected[1]], NaN
            else
                return NaN, Is[selected[1]]
            end
        else
            return NaN, NaN
        end
    end
    
    Icu = getindex.(Ics, 1)
    Icl = getindex.(Ics, 2)
    
    Ic = map(zip(Icu, Icl)) do (u, l)
        if abs(u + l) > 20
            return NaN
        end
        return (u - l) / 2
    end
    
    return Ic, Icu, Icl
end

function compute_normal_resistances(dVdImat::Matrix, Is::Vector, Ic::Vector)
    """Compute normal-state resistances."""
    RNs = map(axes(dVdImat, 1)) do ir
        if !isnan(Ic[ir]) && Ic[ir] != 0
            threshold = abs(Ic[ir]) + 20
        else
            threshold = 20
        end
        iIs = findall(x -> abs(x) >= threshold, Is)
        if !isempty(iIs)
            return mean(dVdImat[ir, iIs])
        else
            return NaN
        end
    end
    return RNs
end

function prepare_data(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120)
    """Load and preprocess experimental data."""
    data = read_data(exp_data)
    @unpack Ts, Is, dVdImat, dVmat, Vmat = data 
    
    Is .*= 1e9
    idx = findall(Is .<= Imax .&& Is .>= -Imax) 
    Is = Is[idx]
    dVdImat = dVdImat[:, idx]
    dVmat = dVmat[:, idx]
    Vmat = Vmat[:, idx]
    
    return Ts, Is, dVdImat, dVmat, Vmat
end

function find_minima_indices(Rf::Vector; min_prom=2)
    """Find minima indices in Rf at I=0."""
    pks = findminima(Rf;)
    pks = peakproms(pks; min = min_prom)
    indeces = pks.indices[2:end] 
    pushfirst!(indeces, 1)
    return indeces
end

function extract_Ic(name::String; basepath = "data/Js")
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng = params
    
    Ic = getindex(findmax(Js; dims = 2),1) |> vec
    return Brng, Ic * 74.4 # in nA
end
# ============================================================================
# Individual figure functions
# ============================================================================

function fig_conductance_heatmap(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120; plotIc::Bool=true)
    """Figure 1: Conductance heatmap with critical currents overlay."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    
    fig = Figure(size = (600, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    
    hmap = heatmap!(ax, Ts, Is, log10.(abs.(dVdImat)); colormap = :thermal, colorrange = (2.9, 3.2), rasterize = 5)
    Colorbar(fig[1, 2], hmap; label = L"\log\frac{dV}{dI}")
    
    plotIc && scatter!(ax, Ts, Ic; color = :black, markersize = 8)
    plotIc && scatter!(ax, Ts, -Ic; color = :black, markersize = 8)
    
    return fig
end

function fig_resistance_I0(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120)
    """Figure 2: dV/dI vs B at I=0 with colored minima."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    Rf = savitzky_golay(R0, 7, 1).y
    indeces = find_minima_indices(Rf)
    
    colors = (:red, :blue, :green, :orange, :purple)
    
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$dV/dI$ ($\Omega$)",
    )
    
    lines!(ax, Ts, R0; color = :black, linewidth = 2, linestyle = :dash)
    lines!(ax, Ts, Rf; color = :black)
    
    for (i, color) in zip(indeces, colors)
        scatter!(ax, Ts[i], Rf[i]; color = color, markersize = 15)
    end
    
    xlims!(ax, minimum(Ts), maximum(Ts))
    
    return fig, indeces
end

function fig_conductance_currents(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120)
    """Figure 3: dV/dI vs I with peak analysis."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    Rf = savitzky_golay(R0, 7, 1).y
    indeces = find_minima_indices(Rf)
    
    colors = (:red, :blue, :green, :orange, :purple)
    
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$dV/dI$ ($\Omega$)",
        ylabel = L"$I$ (nA)",
    )
    
    ylims!(ax, minimum(Is), maximum(Is))
    xlims!(ax, 0, maximum(dVdImat[1, :]) * 1.1)
    
    for (i, color) in zip(indeces, colors)
        lines!(ax, dVdImat[i, :], Is; color = color, linewidth = 2)
        all_indices = get_peak_indices_smoothed(dVdImat[i, :])
        selected = select_top_peaks(all_indices, Is, 2)
        scatter!(ax, dVdImat[i, selected], Is[selected]; color = color, markersize = 10)
    end
    
    # Background gray peaks for all B fields
    for ir in axes(dVdImat, 1)
        all_indices = get_peak_indices_smoothed(dVdImat[ir, :])
        selected = select_top_peaks(all_indices, Is, 2)
        scatter!(ax, dVdImat[ir, selected], Is[selected]; color = :gray, alpha = 0.3, markersize = 5)
    end
    
    return fig
end

function fig_critical_current(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120)
    """Figure 4: Critical current vs B."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c$ (nA)",
    )
    
    scatter!(ax, Ts, Ic; color = :black, markersize = 8)
    xlims!(ax, minimum(Ts), maximum(Ts))
    
    return fig, Ic
end

function fig_normal_resistance(exp_data::Exp_Data, exp_params::Exp_Params, Imax::Real=120, Blims::Tuple=(0.11, 0.4))
    """Figure 5: Normal resistance vs B with excluded regions."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Ic)
    
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_N$ ($\Omega$)",
    )
    
    scatter!(ax, Ts, RNs; color = :black, markersize = 8)
    xlims!(ax, minimum(Ts), maximum(Ts))
    band!(ax, [0, Blims[1]], minimum(RNs), maximum(RNs); color = (:gray, 0.2))
    band!(ax, [Blims[2], maximum(Ts)], minimum(RNs), maximum(RNs); color = (:gray, 0.2))
    
    return fig, RNs
end

function fig_rnorm_and_ic_ah(exp_data::Exp_Data, exp_params::Exp_Params, T::Real=250, Imax::Real=120, Blims::Tuple=(0.11, 0.4))
    """Figure 6: Normalized resistance and Andreev-Huang critical current."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Ic)
    
    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]
    
    Rnorm = R0 ./ RNs
    γs = γAH.(Rnorm)
    IcAH(T_val) = 0.0209 * T_val .* γs  # T in mK
    
    fig = Figure(size = (900, 400), fontsize = 20)
    
    ax1 = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_p / R_N$",
    )
    scatter!(ax1, Ts[ilim:imax], Rnorm[ilim:imax]; color = :black, markersize = 8)
    xlims!(ax1, minimum(Ts), maximum(Ts))
    
    ax2 = Axis(fig[1, 2];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c^{AH}$ (nA)",
    )
    scatter!(ax2, Ts[ilim:imax], IcAH(250)[ilim:imax]; color = :red, markersize = 8, label = L"$T=250\,\mathrm{mK}$")
    scatter!(ax2, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :green, markersize = 8, label = L"$T=%$(T)\,\mathrm{mK}$")
    xlims!(ax2, minimum(Ts), maximum(Ts))
    ylims!(ax2, 0, maximum(filter(!isnan, IcAH(T))) * 1.1)
    axislegend(ax2, position = :rt)
    
    return fig, IcAH
end

function fig_zoom_comparison(exp_data::Exp_Data, exp_params::Exp_Params, T::Real=250, Imax::Real=120, IM::Real=40, Blims::Tuple=(0.11, 0.4))
    """Figure 7: Zoomed heatmap with critical current comparison."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Ic)
    
    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    
    Rnorm = R0 ./ RNs
    γs = γAH.(Rnorm)
    IcAH(T_val) = 0.0209 * T_val .* γs  # T in mK
    
    izoom = findall(Is .<= IM .&& Is .>= -IM) 
    
    fig = Figure(size = (600, 400), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    
    heatmap!(ax, Ts, Is[izoom], log10.(abs.(dVdImat[:, izoom])); colormap = :thermal, colorrange = (2.9, 3.2), rasterize = 5)
    scatter!(ax, Ts[ilim:imax], IcAH(250)[ilim:imax]; color = :red, markersize = 8, label = L"$T=250\,\mathrm{mK}$")
    scatter!(ax, Ts[ilim:imax], -IcAH(250)[ilim:imax]; color = :red, markersize = 8)
    scatter!(ax, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :green, markersize = 8, label = L"$T=%$(T)\,\mathrm{mK}$")
    scatter!(ax, Ts[ilim:imax], -IcAH(T)[ilim:imax]; color = :green, markersize = 8)
    
    axislegend(ax, position = :rt)
    
    fullIc = vcat(Ic[1:ilim-1], IcAH(T)[ilim:imax])
    while length(fullIc) < length(Ts)
        push!(fullIc, NaN)
    end
    return fig, fullIc, Ts
end

# ============================================================================
# Combined figure (original layout)
# ============================================================================

function fig_data(exp_data::Exp_Data, exp_params::Exp_Params, T::Real=250; 
                 colors::Tuple=(:red, :blue, :green, :orange, :purple), 
                 Imax::Real=120, 
                 Blims::Tuple=(0.11, 0.4))
    """Combined figure with all panels."""
    
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    Rf = savitzky_golay(R0, 7, 1).y
    indeces = find_minima_indices(Rf)
    
    Ic, Icu, Icl = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Ic)
    
    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]
    
    Rnorm = R0 ./ RNs
    γs = γAH.(Rnorm)
    IcAH(T_val) = 0.0209 * T_val .* γs  # T in mK
    
    fig = Figure(size = (900, 900), fontsize = 20)
    
    # Panel 1: Conductance heatmap
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    hmap = heatmap!(ax, Ts, Is, log10.(abs.(dVdImat)); colormap = :thermal, colorrange = (2.9, 3.2), rasterize = 5)
    Colorbar(fig[1, 3], hmap; label = L"\log\frac{dV}{dI}")
    hidexdecorations!(ax, ticks = false)
    scatter!(ax, Ts, Ic; color = :black, markersize = 8)
    scatter!(ax, Ts, -Ic; color = :black, markersize = 8)
    
    # Panel 2: R at I=0
    axR = Axis(fig[2, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$dV/dI$ ($\Omega$)",
    )
    lines!(axR, Ts, R0; color = :black, linewidth = 2, linestyle = :dash)
    lines!(axR, Ts, Rf; color = :black)
    xlims!(axR, minimum(Ts), maximum(Ts))
    hidexdecorations!(axR, ticks = false)
    for (i, color) in zip(indeces, colors)
        scatter!(axR, Ts[i], Rf[i]; color = color, markersize = 15)
    end
    Label(fig[2, 1, Top()], L"$I = 0$"; padding = (200, 0, -220, 0))
    
    # Panel 3: Conductance vs I
    axI = Axis(fig[1, 2];
        xlabel = L"$dV/dI$ ($\Omega$)",
        ylabel = L"$I$ (nA)",
    )
    ylims!(axI, minimum(Is), maximum(Is))
    xlims!(axI, 0, maximum(dVdImat[1, :]) * 1.1)
    hideydecorations!(axI, ticks = false)
    for (i, color) in zip(indeces, colors)
        lines!(axI, dVdImat[i, :], Is; color = color, linewidth = 2)
        all_indices = get_peak_indices_smoothed(dVdImat[i, :])
        selected = select_top_peaks(all_indices, Is, 2)
        scatter!(axI, dVdImat[i, selected], Is[selected]; color = color, markersize = 10)
    end
    for ir in axes(dVdImat, 1)
        all_indices = get_peak_indices_smoothed(dVdImat[ir, :])
        selected = select_top_peaks(all_indices, Is, 2)
        scatter!(axI, dVdImat[ir, selected], Is[selected]; color = :gray, alpha = 0.3, markersize = 5)
    end
    
    # Panel 4: Critical current
    axIc = Axis(fig[3, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c$ (nA)",
    )
    scatter!(axIc, Ts, Ic; color = :black, markersize = 8)
    xlims!(axIc, minimum(Ts), maximum(Ts))
    hidexdecorations!(axIc, ticks = false)
    
    # Panel 5: Normal resistance
    axRN = Axis(fig[3, 2];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_N$ ($\Omega$)",
    )
    scatter!(axRN, Ts, RNs; color = :black, markersize = 8)
    xlims!(axRN, minimum(Ts), maximum(Ts))
    hidexdecorations!(axRN, ticks = false)
    band!(axRN, [0, Blims[1]], minimum(RNs), maximum(RNs); color = (:gray, 0.2))
    band!(axRN, [Blims[2], maximum(Ts)], minimum(RNs), maximum(RNs); color = (:gray, 0.2))
    
    # Panel 6: Rnorm
    axRnorm = Axis(fig[4, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$R_p / R_N$",
    )
    scatter!(axRnorm, Ts[ilim:imax], Rnorm[ilim:imax]; color = :black, markersize = 8)
    xlims!(axRnorm, minimum(Ts), maximum(Ts))
    hidexdecorations!(axRnorm, ticks = false)
    
    # Panel 7: IcAH
    axIcAH = Axis(fig[4, 2];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c^{AH}$ (nA)",
    )
    scatter!(axIcAH, Ts[ilim:imax], IcAH(250)[ilim:imax]; color = :red, markersize = 8)
    scatter!(axIcAH, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :green, markersize = 8)
    scatter!(axIc, Ts[ilim:imax], IcAH(250)[ilim:imax]; color = :red, markersize = 8)
    xlims!(axIcAH, minimum(Ts), maximum(Ts))
    ylims!(axIcAH, 0, maximum(filter(!isnan, IcAH(250))) * 1.1)
    
    # Panel 8: Zoom
    axZoom = Axis(fig[5, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )
    IM = 40
    izoom = findall(Is .<= IM .&& Is .>= -IM) 
    heatmap!(axZoom, Ts, Is[izoom], log10.(abs.(dVdImat[:, izoom])); colormap = :thermal, colorrange = (2.9, 3.2), rasterize = 5)
    scatter!(axZoom, Ts[ilim:imax], IcAH(250)[ilim:imax]; color = :red, markersize = 8, label = L"$T = 250\,\mathrm{mK}$")
    scatter!(axZoom, Ts[ilim:imax], -IcAH(250)[ilim:imax]; color = :red, markersize = 8)
    scatter!(axZoom, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :green, markersize = 8, label = L"$T = %$(T)\,\mathrm{mK}$")
    scatter!(axZoom, Ts[ilim:imax], -IcAH(T)[ilim:imax]; color = :green, markersize = 8)
    Legend(fig[5, 2], axZoom; position = (0, 0), framevisible = false)
    
    colsize!(fig.layout, 2, Relative(0.4))
    rowgap!(fig.layout, 1, -40)
    rowgap!(fig.layout, 2, 5)
    rowgap!(fig.layout, 4, 5)
    rowgap!(fig.layout, 4, -40)
    colgap!(fig.layout, 2, 5)
    
    return fig
end

# ============================================================================
# Theory comparison
# ============================================================================
##
function fig_theory(exp_data::Exp_Data, exp_params::Exp_Params, name::String; name_base = "C_mu=25_w=20_tau=2.5_g=17_T=550_TN=0.19", T::Real = 250, Imax::Real=120, Blims = (0.11, 0.4), pref = 0.8, pref_base = 0.8)
    """Figure: Comparison of experimental and theoretical critical currents."""
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, exp_params, Imax)
    
    Ic, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Ic)
    
    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]
    
    i0 = findmin(abs.(Is))[2]
    R0 = dVdImat[:, i0]
    
    Rnorm = R0 ./ RNs
    γs = γAH.(Rnorm)
    IcAH(T_val) = 0.0209 * T_val .* γs  # T in mK


    
    fig = Figure(size = (600, 800), fontsize = 20)
    ax = Axis(fig[1, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I$ (nA)",
    )

    heatmap!(ax, Ts, Is, log10.(abs.(dVdImat)); colormap = :thermal, colorrange = (2.9, 3.2), rasterize = 5)
    
    Brng_base, Ict_base = extract_Ic(name_base)
    Brng, Ict = extract_Ic(name)
    scatter!(ax, Ts[1:ilim], Ic[1:ilim]; color = :black, markersize = 8)
    scatter!(ax, Ts[1:ilim], -Ic[1:ilim]; color = :black, markersize = 8)

    scatter!(ax, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :black, markersize = 8, label = "Experiment")
    scatter!(ax, Ts[ilim:imax], -IcAH(T)[ilim:imax]; color = :black, markersize = 8)

    #lines!(ax, Brng_base, pref_base .* Ict_base; color = (:red, 0.7),  linestyle = :dash, linewidth = 3, label = "Theory (base)")
    #lines!(ax, Brng_base,- pref_base .* Ict_base; color = (:red, 0.7),  linestyle = :dash, linewidth = 3,)
    lines!(ax, Brng, pref .* Ict; color = (:red, 0.7), linewidth = 3, label = "Theory")
    lines!(ax, Brng, -pref .* Ict; color = (:red, 0.7), linewidth = 3)

    axislegend(ax, position = :rt)

    text!(ax, 0.35, -100; text = L"$T=%$(T)\mathrm{mK}$", align = (:center, :center), fontsize = 24, color = :white)

    hidexdecorations!(ax, ticks = false)

    ax = Axis(fig[2, 1];
        xlabel = L"$B$ (T)",
        ylabel = L"$I_c$ (nA)",
    )
    
    scatter!(ax, Ts[1:ilim], Ic[1:ilim]; color = :black, markersize = 8)
    scatter!(ax, Ts[ilim:imax], IcAH(T)[ilim:imax]; color = :black, markersize = 8, label = "Experiment")
    #lines!(ax, Brng_base, pref_base .* Ict_base; color = (:red, 0.7),  linestyle = :dash, linewidth = 3, label = "Theory (base)")
    lines!(ax, Brng, pref .* Ict; color = :red, linewidth = 3, label = "Theory")
    xlims!(ax, minimum(Ts), maximum(Ts))
    ylims!(ax, 0, maximum(filter(!isnan, Ic)) * 1.1)
    
    return fig
end
fig8 = fig_theory(Exp_Data(), Exp_Params(), "fitted_C3"; name_base = "C_mu=25_w=20_tau=2.5_g=17_T=550_TN=0.19", T = 550, pref = 0.8)
fig8
#save("plots/study_figs/fig8_C.pdf", fig8)

##
# ============================================================================
# Usage examples
# ============================================================================

# Combined figure with all panels
fig = fig_data(Exp_Data(), Exp_Params(), 550)

# Individual figures
fig0 = fig_conductance_heatmap(Exp_Data(), Exp_Params(); plotIc=false)
save("plots/study_figs/fig0.pdf", fig0)
fig1 = fig_conductance_heatmap(Exp_Data(), Exp_Params())
save("plots/study_figs/fig1.pdf", fig1)
fig2, indeces = fig_resistance_I0(Exp_Data(), Exp_Params())
fig2
save("plots/study_figs/fig2.pdf", fig2)
fig3 = fig_conductance_currents(Exp_Data(), Exp_Params())
save("plots/study_figs/fig3.pdf", fig3)
fig4, Ic = fig_critical_current(Exp_Data(), Exp_Params())
fig4
save("plots/study_figs/fig4.pdf", fig4)
fig5, RNs = fig_normal_resistance(Exp_Data(), Exp_Params())
fig5
save("plots/study_figs/fig5.pdf", fig5)
fig6, IcAH = fig_rnorm_and_ic_ah(Exp_Data(), Exp_Params(), 500)
fig6
save("plots/study_figs/fig6.pdf", fig6)
fig7, Ic, Ts = fig_zoom_comparison(Exp_Data(), Exp_Params(), 500)
save("plots/figures/Ic_fitting_AH.pdf", fig7)

fig8 = fig_theory(Exp_Data(), Exp_Params(), "fitted_C"; T = 550, pref = 0.85)
save("plots/figures/Ic_fitting_AH_vs_theo.pdf", fig8)
fig8