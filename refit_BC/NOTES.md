# Devices B and C, m_J-resolved fit with the physical Zeeman (started 2026-11-05)

Carlos: "as a possibility for the correction, in a new branch, write a new code (using the m_J
independent Hamiltonian builder) to fit the critical current of devices B and C, with no lobe 2
skewness. If the data is not in this repo, extract it from the pdfs of the figures for a first
estimation (I'll provide the experimental data later) … be as loyal as possible to the original code
and spirit." Then: SOC = 0, g realistic between 0 and 20; the final report reproduces the theory
figure (Fig. 3) for B and C; clusters allowed (esbirro free, atto partly); T_e is a fit parameter
whose only role is relating I_c to R_0. Branch `refit_devices_BC`. Background: `wip/refit_zeeman/`.

## Pipeline

| step | script | env | output |
|---|---|---|---|
| 1 digitize SI figures | `digitize.jl` (+ `plots/kitchen/pdf_digitizer.jl`) | plots | `data_exp/digitized/` |
| 2 experimental targets | `target.jl` (+ `plots/kitchen/digitized_reader.jl`) | plots | `out/target_{A,B,C}.*` |
| 3 grid scan | `scan.jl`, `scan.sbatch` (+ `src/calculations/Ic_single_mode.jl`) | main | `out/<grid>_<tag>.jld2` |
| 4 scoring | `analyze.jl` | plots | `out/<grid>_<mode>_{ranking,per_g}.txt`, `_scores.jld2` |
| 5 final curves | `final.jl` | main | `out/final_<grid>_<mode>_<rank>.jld2` |
| 6 figure | `fig_theo_exp_BC.jl` | plots | `figures/fig_theo_exp_{B,C}.pdf` |

Device geometry and windows: `devices.jl`. Checks: `validate.jl` (device A), `smoke_check.jl`.

## Log

- Digitization: maps and R_0 lines of SI_DEVICE_B (a, c) and SI_DEVICE_C (a, b, c); colour inversion
  through the embedded colourbars (7.7–7.9 Ω levels), tick fits ≤ 9·10⁻⁶ T. B map vs published R_0
  line: rms 16 Ω. C zoom map vs published line: rms 76 Ω, up to 200 Ω in 2L (the line matches the
  full-range map better; caption says it was taken from the zoom). Details `data_exp/digitized/README.md`.
- Lobe edges from the high-bias dV/dI maxima: B 0.080, ≈ 0.238 T → R_LP = 64.5 nm; C 0.060, 0.180,
  0.297, 0.414 T → P = 0.118 T → R_LP = 75 nm (SM text: ~77 nm, ~110 mT). Tilts, ξ_d, d, Δ0 of device A.
- Targets: current offset −3.2 nA removed for B (fig_gates.jl procedure), −0.05 nA for C. T_e from
  I_c^AH = I_sw at the 1L maximum (the paper's procedure; 475 mK on A vs the published 460): B 604 mK,
  C 292 mK. Whole-1L least squares would give 543 mK on A, so it is not the paper's procedure.
- Validation on A (`out/validate.txt`): single-sector τ 0.86242 vs multimode 0.86241; Σ_Z I_c at
  g_1D = 17π within 0.34 % of the published fitted_C3 peak; physical g = 17 lowers 2L ≈ 3×.
- Smoke (es1 2368, 2369): sectors outside |Z| ≤ 9 (B) and ≤ 12 (C) carry ≤ 1.8·10⁻⁴ of I_c at the
  largest μ and R_av; 1.6–1.7 s per (B, Z) point.
- Full scans: B (500 candidates: g 0–20 step 5, w 15–35, μ 10–30, τΓ 2–3.5; T = 600 mK; 60 fields;
  Z −9:9) on esbirro 2370 (8 × 32). C (600 candidates: w 15–40; T = 290 mK; 70 fields to 0.42 T;
  Z −12:12) on atto 262181 (long, 2 × 96, part 1/2) and 262182 (most, 2 × 96, part 2/2).
  262181 was still loading Julia after 15 min (both jobs shared node a4 and loaded 384 workers at
  once); cancelled and its half resubmitted as parts 1/4 (262183, a6) and 3/4 (262184, a2), 96 tasks
  each on otherwise idle nodes. Real cost ≈ 4–6 s per (B, Z) point with full nodes (smoke: 1.6 s on a
  lightly loaded node).
- End-to-end test (`testB`, one candidate: g 15, w 25, μ 25, τΓ 3, local): scan → analyze → final →
  figure all run; this device-A-like candidate already gives rmsR_1L 0.049 and off1 +0.022 T (data +0.017).
- B scan done (es1 2370, 2.09 h on 256 workers). T_e fixed (600 mK): best g 10, w 30 (R_av 45.0 nm),
  μ̃ 15 meV, τΓ 3.0, T_N ≈ 0.20; cost 0.067 (rmsR_1L 0.035, rel0L 0.064), off1 +0.019 T (data +0.017).
  Best per g: cost 0.084 (g 0), 0.074 (5), 0.067 (10), 0.078 (15), 0.076 (20), all with R_av 42.5–45,
  μ̃ 15–20, τΓ 3–3.5 and the right 1L skew: device B does not constrain g. T_e free: cost 0.059 at
  T_e ≈ 1.1 K (R_av 50, μ̃ 20, τΓ 3.5), a marginal gain at an implausible T_e, so T_e stays fixed.
- B2 refinement (81 candidates, half steps around the best) on es1 2371; final curves of the coarse
  best computed locally (`final.jl B fixed 1`).
- B2 done (0.3 h): best g 7.5, w 32.5 (R_av 43.75 nm), μ̃ 15, τΓ 3.25, T_N ≈ 0.20, cost 0.060
  (rmsR_1L 0.031, rel0L 0.058), off1 +0.0179 T (data +0.0173). g 10 equivalent (0.061). Final curves on
  es1 2372 → `figures/fig_theo_exp_B2.pdf`. analyze.jl now uses the target's own T_e in mode "fixed"
  (604 / 292 / 299 mK) instead of the rounded kernel value.
- C scan done (atto 262182–262184, 2.6 h): poor. Best cost 0.36 (g 0, R_av 56, μ̃ 15, τΓ 3.5 = grid
  edge, T_N ≈ 0.11): rmsR 1L 0.19, 2L 0.13, 3L 0.012, rel0L 0.057. The data's 1L is flat (I_c^AH 23–26
  nA, R_0/R_N 0.16–0.22 across the lobe); every candidate's I_c rises from ≈ 10 nA at the left edge to
  ≈ 45 nA at 0.14 T. 2L theory dip centred (off2 −0.001 at g 0, +0.006 at g 20) vs data −0.011. T_e
  free: 0.34. Alternative target `Cline` (published R_0 line / R_N of map a, T_e 299 mK, 1L right-skewed
  like A and B): 0.32, rmsR_1L still 0.18. Cost decreases towards τΓ ↑ and μ̃ ↓ (edges of the grid).
- C2 (144 candidates: τΓ 4, 5, 6; μ̃ 2.5–10; w 25–40; g 0, 10, 20; Z −10:10): atto 262186–262189 (a1,
  a2, a6, a4) and es1 2373 (es1–3). Two esbirro submissions whose batch host was not es1 failed with
  JobLaunchFailure (2374 es4–6, 2375 es4–6, 2376 es5,7,8); that part went to atto a4.
  Result (all 144): best 0.356 (g 10, w 35 = R_av 53.5, μ̃ 10, τΓ 4, T_N ≈ 0.12): rmsR 1L 0.15, 2L 0.12,
  3L 0.044, rel0L 0.072. No real gain over the coarse grid. Node a1 ran ≈ 3× slower than a2/a6.
- C3 (ξ_d of both leads 30, 45, 60 nm; g 10; w 25–35; μ̃ 10–20; τΓ 2.5–3.5; 81 candidates): atto
  262191–262193 (a2, a5, a6), es1 2377 (es1–3 with es1 as batch host: launches fine). Best ξ_d 60, cost
  0.37; ξ_d 45: 0.42; 30 worse. 1L rms 0.15–0.18 at every ξ_d, 2L rms 0.19–0.22. Shorter ξ_d (larger edge
  gaps) does not supply the missing supercurrent near the lobe edges: hypothesis refuted.
- Final curves of the C2 best (atto 262194, a5) → `figures/fig_theo_exp_C.pdf` (renamed from _C2);
  B2 best → `figures/fig_theo_exp_B.pdf` (renamed from _B2). Report: `REPORT.md` / `REPORT.pdf`.
