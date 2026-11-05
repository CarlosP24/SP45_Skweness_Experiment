# Fit target for the Zeeman-corrected refit: the measured R_0/R_N(B) and I_c^AH(B) of the main
# device, extracted exactly as Fig. 3 of the preprint (plots/plotters/theo_exp_core.jl,
# `process_data`, experimental part only): run 250mK_DEV13_275_16 (Exp_Data() default),
# R_s = 5490 Ω, β = 0.95, Savitzky-Golay window 50, B ∈ (0.08, 0.4) T, T_e = 460 mK.
# Run from the SP45 root:  julia --project=plots wip/refit_zeeman/target.jl
using Parameters, JLD2, DelimitedFiles, Statistics, SavitzkyGolay, Roots, Bessels, Peaks
using CairoMakie
include(joinpath(pwd(), "plots", "kitchen", "data_reader.jl"))
include(joinpath(pwd(), "plots", "kitchen", "Ic_reader.jl"))

const Te = 460.0          # mK, the published T_e (from the resolved 1L peak; experimental, kept)
const Blims = (0.08, 0.4)

exp_data = Exp_Data()
Ts, Is, dVdImat, _, _ = prepare_data(exp_data, 120)
Isw, _, _ = extract_peak_currents(dVdImat, Is, Ts)      # dV/dI peak (switching) current, nA
RNs = compute_normal_resistances(dVdImat, Is, Isw)
i0 = findmin(abs.(Is))[2]
r0 = vec(dVdImat[:, i0] ./ RNs)
γ = γAH.(abs.(r0))
IcAH = 0.0209 .* Te .* γ                                   # nA (as published)

iB = findall(B -> Blims[1] < B <= Blims[2], Ts)
B = Ts[iB]
jldsave("wip/refit_zeeman/target.jld2"; B, r0 = r0[iB], IcAH = IcAH[iB], Isw = Isw[iB],
    Ball = Ts, Iswall = Isw, Te, run = exp_data.lockin)
open("wip/refit_zeeman/target.csv", "w") do io
    println(io, "B_T,R0_over_RN,IcAH_nA,Isw_nA")
    for i in eachindex(B)
        println(io, join((B[i], r0[iB][i], IcAH[iB][i], Isw[iB][i]), ","))
    end
end

fig = Figure(size = (700, 450))
ax1 = Axis(fig[1, 1]; ylabel = "I (nA)", title = "fit target, run 275_16, T_e = 460 mK")
scatter!(ax1, Ts, Isw; color = :gray, markersize = 5, label = "dV/dI peak (I_sw)")
scatter!(ax1, B, IcAH[iB]; color = :green, markersize = 6, label = "I_c^AH")
axislegend(ax1; position = :rt)
ax2 = Axis(fig[2, 1]; xlabel = "B (T)", ylabel = "R0/RN")
scatter!(ax2, B, r0[iB]; color = :black, markersize = 6)
linkxaxes!(ax1, ax2)
save("wip/refit_zeeman/target.png", fig)
println("target: $(length(B)) field points, IcAH max $(round(maximum(filter(isfinite, IcAH[iB])); digits = 1)) nA")
