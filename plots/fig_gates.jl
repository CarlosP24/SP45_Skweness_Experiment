include("plotters/theo_exp_core.jl")

# Referee: gate-voltage independence. Measurements at lower V_g than the main-text one (V_g = 80 V).
# The 43_43 run only has the Y quadrature of the lock-in. bad_B: corrupted field sweeps, masked.
const gate_runs = [
    (Vg = 20, run = "A198_DEV19_54_54", quad = "X", bad_B = Float64[]),
    (Vg = 25, run = "A198_DEV19_51_51", quad = "X", bad_B = Float64[]),
    (Vg = 30, run = "A198_DEV19_43_43", quad = "Y", bad_B = Float64[]),
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

function fig_gates(name::String; R = 64, Imax = 120, Bmax = 0.4, pref = 0.8, rrange = (0.5, 1.6), δI_ref = 0.2167)
    Φ0 = 2067.833 #T⋅nm^2
    Bns = [1, 2].* Φ0 / (π * R^2)
    Bhalfs = [0.5, 1.5].* Φ0 / (π * R^2)

    Brng, Ict = ics("fitted_$(name)")
    Ict .*= pref
    Ict0 = Ict[findmin(abs.(Brng))[2]]

    fig = Figure(size = (1100, 560), fontsize = 20)

    axs = Axis[]
    for (i, (; Vg, run, quad, bad_B)) in enumerate(gate_runs)
        row, col = fldmod1(i, 2)

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
        ibad = findall(B -> any(isapprox(B, b; atol = 1e-6) for b in bad_B), Ts)
        dVdImat[ibad, :] .= NaN

        # Theory scaled to the measured zero-field switching current (smaller transparency at lower V_g)
        scale = switching_current(Ts, Is, dVdImat) / Ict0

        # Same colour range relative to R_N as in the main-text figure ((0.5, 1.5) kΩ for R_N ≈ 1 kΩ)
        RN = median(filter(isfinite, dVdImat[Ts .> 0.1, abs.(Is) .> 80]))
        colorrange = round.(rrange .* RN; digits = -2)

        ax = Axis(fig[row, 2col - 1];
            xlabel = L"$B$ (T)",
            ylabel = L"$I$ (nA)",
        )
        lims, hmap = panel_dVdI(ax, dVdImat, Ts, Is, Ts, Brng, Ict .* scale; colormap = :seismic, lowclip = cgrad(:seismic)[1], colorrange, highclip = cgrad(:seismic)[end], legend = false, show_AH = false)
        ax.ylabelpadding = -18

        vlines!(ax, Bns; color = :black, linestyle = :dash)
        vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

        for n in 0:2
            text!(ax, (n + 0.3 - (n == 2) * 0.1) * Φ0 / (π * R^2), -110; text = "$n"*"L", align = (:center, :center), color = :white)
        end

        text!(ax, 0.395, 110; text = L"V_g = %$(Vg)\,\mathrm{V}", align = (:right, :top), color = :white)
        text!(ax, 0.395, 80; text = L"I_c^\text{th} \times %$(round(scale; digits = 2))", align = (:right, :top), color = :white, fontsize = 16)

        Colorbar(fig[row, 2col], hmap; label = L"$dV/dI$ (k$\Omega$)", labelpadding = -20, ticks = ([first(colorrange), last(colorrange)], [first(colorrange), last(colorrange)] ./ 1e3 .|> string))

        row == 1 && hidexdecorations!(ax, ticks = false, grid = false)
        col == 2 && hideydecorations!(ax, ticks = false, grid = false)
        push!(axs, ax)
    end
    linkaxes!(axs...)

    colgap!(fig.layout, 1, 5)
    colgap!(fig.layout, 3, 5)
    rowgap!(fig.layout, 5)

    style = (font = "CMU Serif Bold", fontsize = 20)
    for (i, l) in enumerate(["a", "b", "c", "d"])
        row, col = fldmod1(i, 2)
        Label(fig[row, 2col - 1, TopLeft()], l, padding = (0, 0, -12, 0), halign = :left; style...)
    end

    return fig
end
