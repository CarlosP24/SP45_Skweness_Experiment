# Step 2 of the devices-B/C fit: the experimental targets, extracted as Fig. 3 of the preprint does for
# the main device (plots/plotters/theo_exp_core.jl, `process_data`, experimental part):
#   dV/dI peak current I_sw (extract_peak_currents), R_N (compute_normal_resistances, mean dV/dI over
#   I_sw + 20 nA ≤ |I| ≤ Imax), R_0/R_N at zero bias, I_c^AH = 0.0209 T_e γ_AH(R_0/R_N).
# Two additions, both needed for the SI maps and both already used elsewhere in this repo:
#   * the current-bias offset is removed first (fig_gates.jl, `current_offset`: centre of the 0L switching
#     peaks); device B's map has its zero-bias dip and switching peaks centred at ≈ −3 nA;
#   * T_e is *fitted*, as the paper describes ("fitting I_c^AH to the clearly resolved I_c peak in 1L"):
#     I_c^AH = 0.0209 T_e γ_AH matched to I_sw at the 1L maximum of γ_AH (least squares over that field and
#     its two neighbours). On the main device's run 275_16 this gives 467 mK, the published 460 mK.
#     (A least-squares fit over the whole 1L gives 543 mK there: away from the maximum the dV/dI peak
#     no longer tracks I_c, it stays at 25–30 nA while R_0/R_N → 1.)
#   * T_e is a fit parameter whose only role is this I_c ↔ R_0/R_N relation (Carlos, 2026-11-05), so
#     each device takes its own fitted value (`Te_used = nothing` in devices.jl); device A, a reference
#     here, keeps the published 460 mK.
# Run from the SP45 root:  julia --project=plots refit_BC/target.jl
# Output: refit_BC/out/target_{A,B,C}.jld2, refit_BC/out/target_{B,C}.csv, refit_BC/figures/target.png
using Parameters, JLD2, DelimitedFiles, Statistics, SavitzkyGolay, Roots, Bessels, Peaks, Printf
using CairoMakie
include(joinpath(pwd(), "plots", "kitchen", "data_reader.jl"))
include(joinpath(pwd(), "plots", "kitchen", "Ic_reader.jl"))
include(joinpath(pwd(), "plots", "kitchen", "digitized_reader.jl"))
include(joinpath(pwd(), "plots", "fig_gates.jl"))          # current_offset, switching_peaks
include(joinpath(pwd(), "refit_BC", "devices.jl"))          # lobe edges and windows of A, B, C

const OUTDIR = joinpath("refit_BC", "out")
mkpath(OUTDIR); mkpath(joinpath("refit_BC", "figures"))

"""
    te_fit(B, Isw, γ, window)

T_e (mK) matching I_c^AH = 0.0209 T_e γ to the dV/dI peak current at the 1L maximum of γ: least squares
of Σ (0.0209 T_e γ − I_sw)² over that field and its two neighbours (I_c^AH is linear in T_e at fixed γ).
Returns T_e, the field of the maximum and R_0/R_N there is not returned; the relative rms of the match over
the whole window is returned as a diagnostic.
"""
function te_fit(B, Isw, γ, window)
    m = [window[1] <= b <= window[2] && isfinite(i) && isfinite(g) for (b, i, g) in zip(B, Isw, γ)]
    count(m) < 3 && error("too few resolved 1L peaks for the T_e fit")
    idx = findall(m)
    k = idx[argmax(γ[idx])]
    top = [i for i in k-1:k+1 if i in idx]
    T = sum(γ[top] .* Isw[top]) / (0.0209 * sum(γ[top] .^ 2))
    rel = sqrt(mean((0.0209 * T .* γ[m] .- Isw[m]) .^ 2)) / maximum(Isw[m])
    return T, rel, B[k]
end

function process(dev)
    D = DEVICES[dev]
    Ts, Is, dVdImat, _, _ = prepare_data(D.exp, D.Imax + 30)
    iB = findall(Ts .<= D.Bmax_data)
    Ts, dVdImat = Ts[iB], dVdImat[iB, :]
    Ioff = current_offset(Ts, Is, dVdImat; Bmax = D.edges[1] - 0.02)
    if D.remove_offset
        Is = Is .- Ioff
    end
    iI = findall(abs.(Is) .<= D.Imax)
    Is, dVdImat = Is[iI], dVdImat[:, iI]

    Isw, Icu, Icl = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Isw)
    i0 = findmin(abs.(Is))[2]
    r0 = vec(dVdImat[:, i0] ./ RNs)
    if get(D, :r0_line, false)
        # zero-bias resistance from the published vector line instead of the map row
        Bl, R0l = read_R0_line(D.exp)
        r0 = [lin_resample(Bl, R0l, [b])[1] for b in Ts] ./ RNs
    end
    γ = γAH.(abs.(r0))
    Te_fit, rel, Bpk = te_fit(Ts, Isw, γ, D.te_window)
    Te = D.Te_used === nothing ? Te_fit : D.Te_used
    IcAH = 0.0209 .* Te .* γ
    return (; dev, Ts, Is, dVdImat, Isw, Icu, Icl, RNs, r0, γ, IcAH, Te, Te_fit, te_rel = rel, Bpk,
              Ioff, offset_removed = D.remove_offset, I0 = Is[i0])
end

const DEVS = ("A", "B", "C", "Cline")
results = Dict(dev => process(dev) for dev in DEVS)
for dev in DEVS
    r = results[dev]
    k = argmin(abs.(r.Ts .- r.Bpk))
    @printf("device %s: current offset %+.2f nA (%s), zero-bias row at I = %+.3f nA; 1L maximum at %.3f T (R0/RN = %.3f, I_sw = %.1f nA): T_e fit %.0f mK, used %.0f mK\n",
        dev, r.Ioff, r.offset_removed ? "removed" : "kept, as published", r.I0, r.Bpk, r.r0[k], r.Isw[k], r.Te_fit, r.Te)
end
@printf("published T_e of the main device: 460 mK (fit here: %.0f mK)\n", results["A"].Te_fit)

for dev in DEVS
    r = results[dev]
    jldsave(joinpath(OUTDIR, "target_$dev.jld2"); B = r.Ts, Is = r.Is, dVdImat = r.dVdImat, Isw = r.Isw,
        RNs = r.RNs, r0 = r.r0, IcAH = r.IcAH, Te = r.Te, Te_fit = r.Te_fit, Ioff = r.Ioff,
        offset_removed = r.offset_removed, edges = DEVICES[dev].edges, te_window = DEVICES[dev].te_window)
    dev == "A" && continue
    open(joinpath(OUTDIR, "target_$dev.csv"), "w") do io
        println(io, "B_T,R0_over_RN,IcAH_nA,Isw_nA,RN_Ohm")
        for i in eachindex(r.Ts)
            println(io, join((r.Ts[i], r.r0[i], r.IcAH[i], r.Isw[i], r.RNs[i]), ","))
        end
    end
end

fig = Figure(size = (2000, 650), fontsize = 14)
for (k, dev) in enumerate(DEVS)
    r = results[dev]; D = DEVICES[dev]
    ax = Axis(fig[1, k]; ylabel = "I (nA)", title = @sprintf("device %s: I_sw (grey), I_c^AH (green), T_e = %.0f mK (fit %.0f)", dev, r.Te, r.Te_fit))
    scatter!(ax, r.Ts, r.Isw; color = :gray, markersize = 5)
    m = r.Ts .>= D.edges[1]
    scatter!(ax, r.Ts[m], r.IcAH[m]; color = :limegreen, markersize = 6, strokewidth = 0.5)
    vspan!(ax, [D.te_window[1]], [D.te_window[2]]; color = (:orange, 0.15))
    vlines!(ax, D.edges; color = :gray, linestyle = :dash)
    ylims!(ax, 0, 1.15 * maximum(filter(isfinite, r.Isw)))
    axR = Axis(fig[2, k]; xlabel = "B (T)", ylabel = "R0/RN")
    scatter!(axR, r.Ts[m], r.r0[m]; color = :black, markersize = 5)
    vlines!(axR, D.edges; color = :gray, linestyle = :dash)
    ylims!(axR, 0, 1.1)
    linkxaxes!(ax, axR)
end
save(joinpath("refit_BC", "figures", "target.png"), fig)
println("saved refit_BC/figures/target.png")
