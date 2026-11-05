# Devices of the devices-B/C fit, shared by the experimental side (plots environment) and the theory side
# (main environment). The experiment block is only defined when the kitchen readers are loaded.
#
# Shell geometry. Device A (main text) has every lead parameter from Joule spectroscopy (SM Table S1):
# R_LP = 64 nm, d = 8 nm, ξ_d = 55 / 92 nm, θ = 8° / 5°, Δ0 = 0.204 meV (the fitted_C3 wires).
# Devices B and C have no Joule characterisation, so, as a first estimate:
#   * R_LP from the Little–Parks period measured on the digitized maps (lobe edges = maxima of the
#     high-bias dV/dI, where the shell gap is smallest), with the tilts of device A:
#       B: edges 0.080 and ≈ 0.238 T  → P ≈ 0.159 T → R_LP √cosθ = 64.3 nm → R_LP = 64.5 nm
#       C: edges 0.060, 0.180, 0.297, 0.414 T → P = 0.118 T → R_LP √cosθ = 74.7 nm → R_LP = 75 nm
#     (the SM quotes "R_LP ~ 77 nm, ~110 mT" for C; the measured period is 118 mT);
#   * d, ξ_d, θ and Δ0 of device A. Device B was fabricated on the same chip and measured in the same
#     cool-down as A; C is a different chip, so its ξ_d, θ are the least certain inputs of this fit.
# R is the core radius, R_LP = R + d/2.

const Φ0_Tnm2 = 2067.833848

const GEOMETRY = Dict(
    "A" => (R = 60.0, d = 8.0, ξdL = 55.0, ξdR = 92.0, θL = 8.0, θR = 5.0, Δ0 = 0.204),
    "B" => (R = 60.5, d = 8.0, ξdL = 55.0, ξdR = 92.0, θL = 8.0, θR = 5.0, Δ0 = 0.204),
    "C" => (R = 71.0, d = 8.0, ξdL = 55.0, ξdR = 92.0, θL = 8.0, θR = 5.0, Δ0 = 0.204),
)

"field of flux Φ (in Φ0) through the shell of `dev` (left-lead tilt), T"
B_of_Φ(dev, Φ) = (g = GEOMETRY[dev]; Φ * Φ0_Tnm2 / (π * (g.R + g.d / 2)^2 * cosd(g.θL)))

# Measured lobe edges (T) and fit windows. Lobe n spans (edges[n], edges[n+1]); 0L is B < edges[1].
const EDGES = Dict(
    "A" => [B_of_Φ("A", 0.5), B_of_Φ("A", 1.5), B_of_Φ("A", 2.5)],      # as in Fig. 3
    "B" => [0.080, 0.238],
    "C" => [0.060, 0.180, 0.297, 0.414],
)

# T_e is a fit parameter of the Ambegaokar–Halperin relation only (it links I_c to R_0/R_N), fitted per
# device as in the paper; `Te_used = nothing` takes that fit. Device A (reference only) keeps the
# published 460 mK.
if @isdefined(Digitized_Data)
    const DEVICES = Dict(
        # main device, exactly as Fig. 3 (no offset correction there)
        "A" => (exp = Exp_Data(), Imax = 120, Bmax_data = 0.4, edges = EDGES["A"], Te_used = 460.0,
                te_window = (EDGES["A"][1] + 0.005, EDGES["A"][2] - 0.005), remove_offset = false),
        "B" => (exp = Digitized_Data(; map = joinpath("data_exp", "digitized", "devB_dVdI_map.dat"),
                                      line = joinpath("data_exp", "digitized", "devB_R0_line.dat"), Rs = 5490),
                Imax = 120, Bmax_data = 0.24, edges = EDGES["B"], Te_used = nothing,
                te_window = (EDGES["B"][1] + 0.005, EDGES["B"][2] - 0.005), remove_offset = true),
        # C: the low-current zoom (panel b) is the analogue of the main device's run 275_16; it spans
        # ±99 nA, so R_N is averaged up to 95 nA instead of 120 nA
        "C" => (exp = Digitized_Data(; map = joinpath("data_exp", "digitized", "devC_dVdI_zoom.dat"),
                                      line = joinpath("data_exp", "digitized", "devC_R0_line.dat"), Rs = 5350),
                Imax = 95, Bmax_data = 0.48, edges = EDGES["C"], Te_used = nothing,
                te_window = (EDGES["C"][1] + 0.005, EDGES["C"][2] - 0.005), remove_offset = true),
        # C, alternative target: the published R_0 line of panel (c) over R_N of the full-range map (a),
        # which that line matches better than the zoom (digitization_check.png); I_sw from map (a)
        "Cline" => (exp = Digitized_Data(; map = joinpath("data_exp", "digitized", "devC_dVdI_full.dat"),
                                          line = joinpath("data_exp", "digitized", "devC_R0_line.dat"), Rs = 5350),
                    Imax = 120, Bmax_data = 0.48, edges = EDGES["C"], Te_used = nothing, r0_line = true,
                    te_window = (EDGES["C"][1] + 0.005, EDGES["C"][2] - 0.005), remove_offset = true),
    )
end
