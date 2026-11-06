include("plotters/theo_exp_core.jl")

# Referee: gate-voltage independence. Measurements at lower V_g than the main-text one (V_g = 80 V).
# The 43_43 run only has the Y quadrature of the lock-in. bad_B: corrupted field sweeps, masked.
const gate_runs = [
    (Vg = 20, run = "A198_DEV19_54_54", quad = "X", bad_B = Float64[]),
    (Vg = 25, run = "A198_DEV19_51_51", quad = "X", bad_B = Float64[]),
    # (Vg = 30, run = "A198_DEV19_43_43", quad = "Y", bad_B = Float64[]),  # dropped from the figure
    (Vg = 35, run = "A198_DEV19_17_17", quad = "X", bad_B = [0.124]),
]

gate_exp_data(run, quad) = Exp_Data(;
    keysight = "data_exp/$(run)_dmm_keysight_signal.dat",
    lockin = "data_exp/$(run)_mf1_$(quad).dat"
)

function current_step(exp_data::Exp_Data)
    """Current step of the bias sweep (nA)."""
    Is = sort(unique(DelimitedFiles.readdlm(exp_data.keysight)[:, 2]))
    return median(diff(Is)) * 1e9
end

function anchored_dI(exp_data::Exp_Data; Rs = Exp_Params().Rs, Ithr = 2e-9)
    """Effective lock-in excitation fixed by R0 = 0 at the lowest field (supercurrent branch).
    Used when the lock-in quadrature is not calibrated (e.g. Y instead of X)."""
    lockin = DelimitedFiles.readdlm(exp_data.lockin)
    B0 = minimum(lockin[:, 1])
    idx = findall(r -> lockin[r, 1] == B0 && abs(lockin[r, 2]) < Ithr, axes(lockin, 1))
    return mean(abs.(lockin[idx, 4])) / Rs
end

function gate_exp_params(exp_data::Exp_Data, quad::String)
    # These runs use half the lock-in excitation of the main-text ones (2.5 nA instead of 5 nA)
    quad == "X" && return Exp_Params(; dI = 2.5e-9 * Exp_Params().β)
    return Exp_Params(; dI = anchored_dI(exp_data))
end

function switching_peaks(Ts, Is, dVdImat; Bmax = 0.012, Iwin = 110)
    """Positive and negative switching peaks of dV/dI (nA) for each field B ≤ Bmax."""
    ip = findall(0 .< Is .< Iwin)
    im = findall(-Iwin .< Is .< 0)
    rows = findall(Ts .<= Bmax)
    Ips = [Is[ip][argmax(dVdImat[r, ip])] for r in rows]
    Ims = [Is[im][argmax(dVdImat[r, im])] for r in rows]
    return Ips, Ims
end

function current_offset(Ts, Is, dVdImat; Bmax = 0.05)
    """Current-bias offset (nA): centre of the switching peaks across the 0L lobe."""
    Ips, Ims = switching_peaks(Ts, Is, dVdImat; Bmax)
    return median((Ips .+ Ims) ./ 2)
end

function switching_current(Ts, Is, dVdImat; Bmax = 0.012)
    """Zero-field switching current (nA), symmetrised."""
    Ips, Ims = switching_peaks(Ts, Is, dVdImat; Bmax)
    return (median(Ips) - median(Ims)) / 2
end

function zero_bias_resistance(Ts, Is, dVdImat; Blims = (0.08, 0.4))
    """R_0/R_N at zero bias for Blims[1] ≤ B ≤ Blims[2], with R_N as in fig_theo_exp (process_data)."""
    Isw, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Isw)
    i0 = findmin(abs.(Is))[2]
    iB = findall(Blims[1] .<= Ts .<= Blims[2])
    return Ts[iB], dVdImat[iB, i0] ./ RNs[iB]
end

function fig_gates(name::String; R = 64, Imax = 120, Bmax = 0.4, pref = 0.8, rrange = (0.5, 1.6), δI_ref = 0.2167, T = 460, Blims = (0.08, 0.4))
    Φ0 = 2067.833 #T⋅nm^2
    Bns = [1, 2].* Φ0 / (π * R^2)
    Bhalfs = [0.5, 1.5].* Φ0 / (π * R^2)

    Brng, Ict = ics("fitted_$(name)")
    Ict .*= pref
    Ict0 = Ict[findmin(abs.(Brng))[2]]

    jmax = findmin(abs.(Brng .- Blims[2]))[2]
    Brng_Rs = Brng[1:jmax]

    fig = Figure(size = (1240, 560), fontsize = 20)  # page width ≈ 2.07 × the 600-px column figures

    axs, axRs = Axis[], Axis[]
    for (col, (; Vg, run, quad, bad_B)) in enumerate(gate_runs)
        exp_data = gate_exp_data(run, quad)
        exp_params = gate_exp_params(exp_data, quad)

        # Keep the smoothing window of the main-text data (δI_ref nA step) in units of current
        window_size = max(5, round(Int, 50 * δI_ref / current_step(exp_data)))
        Ts, Is, dVdImat, _, _ = prepare_data(exp_data, Imax + 30; exp_params, window_size)

        iB = findall(Ts .<= Bmax)
        Ts, dVdImat = Ts[iB], dVdImat[iB, :]

        # Remove the current-bias offset and crop symmetrically around zero bias
        Ioff = current_offset(Ts, Is, dVdImat)
        Is = Is .- Ioff
        iI = findall(abs.(Is) .<= Imax)
        Is, dVdImat = Is[iI], dVdImat[:, iI]
        @info "V_g = $Vg V: current offset $(round(Ioff; digits = 1)) nA"

        # Zero-bias cut before masking, so the peak finder never sees NaN rows
        Trng_Rs, Rs = zero_bias_resistance(Ts, Is, dVdImat; Blims)
        Rs[[any(isapprox(B, b; atol = 1e-6) for b in bad_B) for B in Trng_Rs]] .= NaN

        ibad = findall(B -> any(isapprox(B, b; atol = 1e-6) for b in bad_B), Ts)
        dVdImat[ibad, :] .= NaN

        # Theory scaled to the measured zero-field switching current (smaller transparency at lower V_g)
        scale = switching_current(Ts, Is, dVdImat) / Ict0
        @info "V_g = $Vg V: theory scale $(round(scale; digits = 3)) (× pref = $(round(scale * pref; digits = 3)))"

        # Theoretical R_0/R_N as in fig_theo_exp (Ambegaokar-Halperin at T), from the scaled I_c^th
        Rs_theo = @. besseli0(scale * Ict[1:jmax] / (0.0209 * T) / 2)^(-2)

        # Same colour range relative to R_N as in the main-text figure ((0.5, 1.5) kΩ for R_N ≈ 1 kΩ)
        RN = median(filter(isfinite, dVdImat[Ts .> 0.1, abs.(Is) .> 80]))
        colorrange = round.(rrange .* RN; digits = -2)

        ax = Axis(fig[1, 2col - 1];
            ylabel = L"$I$ (nA)",
            title = L"V_g = %$(Vg)\,\mathrm{V}",
        )
        lims, hmap = panel_dVdI(ax, dVdImat, Ts, Is, Ts, Brng, Ict .* scale; colormap = :seismic, lowclip = cgrad(:seismic)[1], colorrange, highclip = cgrad(:seismic)[end], legend = false, show_AH = false)
        ax.ylabelpadding = -18

        vlines!(ax, Bns; color = :black, linestyle = :dash)
        vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

        for n in 0:2
            text!(ax, (n + 0.3 - (n == 2) * 0.1) * Φ0 / (π * R^2), -110; text = "$n"*"L", align = (:center, :center), color = :black)
        end


        Colorbar(fig[1, 2col], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20, ticks = ([first(colorrange), last(colorrange)], [first(colorrange), last(colorrange)] ./ 1e3 .|> string))

        col == 1 && axislegend(ax, [[LineElement(color = :black, linewidth = 4), LineElement(color = :orange, linewidth = 3)]], ["Theory"];
            position = :rt, patchsize = (15, 12), padding = (6, 6, 2, 2), margin = (0, 4, 0, 4))

        hidexdecorations!(ax, ticks = false, grid = false)
        col > 1 && hideydecorations!(ax, ticks = false, grid = false)
        push!(axs, ax)

        axR = Axis(fig[2, 2col - 1];
            xlabel = L"$B$ (T)",
            ylabel = L"$R_0/R_N$",
        )
        panel_resistance(axR, Trng_Rs, Rs, Brng_Rs, Rs_theo; legend = false)
        # Compact legend in the empty corner between the 1L-2L data and the lobe labels
        col == 1 && axislegend(axR; position = (1, 0.42), patchsize = (15, 12), padding = (6, 6, 2, 2), rowgap = 0, margin = (0, 4, 0, 0))
        xlims!(axR, lims)

        vlines!(axR, Bns; color = :black, linestyle = :dash)
        vlines!(axR, Bhalfs; color = :gray, linestyle = :dash)

        for n in 0:2
            text!(axR, (n + 0.4 - (n == 2) * 0.1) * Φ0 / (π * R^2), 0.1; text = "$n"*"L", align = (:center, :center), color = :black)
        end

        col > 1 && hideydecorations!(axR, ticks = false, grid = false)
        push!(axRs, axR)
    end
    linkaxes!(axs...)
    linkaxes!(axRs...)
    for (ax, axR) in zip(axs, axRs)
        linkxaxes!(ax, axR)
    end

    rowsize!(fig.layout, 2, Relative(0.35))
    for c in 1:2:5
        colgap!(fig.layout, c, 5)
    end
    rowgap!(fig.layout, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    for (i, l) in enumerate(["a", "b", "c", "d", "e", "f"])
        row, col = fldmod1(i, 3)
        # Columns 2-3 have no tick labels, so push their labels off the axes
        Label(fig[row, 2col - 1, TopLeft()], l, padding = (0, col > 1 ? 25 : 0, -12, 0), halign = :left; style...)
    end

    return fig
end
