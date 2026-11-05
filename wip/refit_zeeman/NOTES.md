# Zeeman-corrected refit of the Paya:26c junction (started 2026-11-05)

Carlos: "re-do the fitting in SP45, using the correct Zeeman and keeping g between 15 and 20";
cluster use approved. Nothing in `src/` has been changed; everything lives in this folder.

## Why

FullShell `build_cyl_mm` (multimode, used for every SP45 wire) writes the Zeeman term as
`0.5 g μBΦ0 πoΦ0 B` = π × ½gμ_B B. `build_cyl` (single mode) gives the physical
½gμ_BΦ0·(πBR²/Φ0)/(πR²) = ½gμ_B B. The published fit (`fitted_C3`, g = 17) therefore has a
physical g = 17π ≈ 53.4. Details: `../../../PT_Thesis/PT5_Josephson/HANDOFF.md`, Q1.

## Method

- **Sectors instead of multimode.** Both leads have R = 60, d = 8 nm (same fluxoid number) and the
  junction is cylindrically symmetric, so m_J is conserved: Σ_Z of single-mode `build_cyl` currents is
  the multimode current, with the physical Zeeman. `ic_single_mode.jl` (`sm_wire`, `js_serial`,
  `ic_single_mode`), built on `src/models/cpr.jl` and `src/calculations/CPR.jl`.
- **Validation** (`validate.jl` → `validate.txt`): at the published convention (g_1D = 17π) the
  sector sum reproduces `data/Js/fitted_C3.jld2` within 0.5 % of the peak at Φ/Φ0 = 0.25–2.0.
  With the physical g = 17 the ratio to the published I_c is 1.01, 0.62, 1.00, 1.14, 0.31, 0.39
  at Φ/Φ0 = 0.25, 0.75, 1.0, 1.25, 1.75, 2.0: lobe 2 drops about threefold. Timing: ≈ 1.3 s per
  (Φ, Z) point per core.
- **Target** (`target.jl` → `target.{jld2,csv,png}`): R_0/R_N and I_c^AH of run `250mK_DEV13_275_16`
  (the one Fig. 3 uses, `Exp_Data()` default), published processing, T_e = 460 mK; plus the dV/dI
  peak current I_sw in lobe 0. Reproduces Fig. 3(c, d).
- **τ**: calibrated per (μ, w, τΓ) to T_N = 0.24 with the hero's multimode normal-state conductance
  at B = 0 (no Zeeman there, so the bug can't enter).
- **Cost** (`analyze.jl`): best overall scale s of I_c (transparency rescaling, T_N ≈ 0.24 s, as
  ch5-plan D11), then rms[R_th − R_0/R_N] over 0.085–0.38 T + 0.5 · rms[(I_c − I_sw)/I_sw(0)] over
  lobe 0, with R_th from the published AH conversion (74.4 nA per code unit, T_e = 460 mK).
- **Baseline**: the published fit scores cost 0.174 (rms R 0.097, lobe-0 0.154) with s = 0.84
  (the paper used 0.8).

## Runs

| Name | Grid | Where | Status |
|---|---|---|---|
| smoke | g 17; w 25; μ 20; τΓ 2.5, 3.0; 4 fields | local, 4 workers, 2.2 min | ok; |Z| ≥ 8 sectors carry ≤ 2·10⁻⁴ of I_c |
| scan1 | g 15, 17.5, 20; w 15–35 (step 5); μ 10–30 (step 5); τΓ 2–3.5 (step 0.5); 76 fields 0.005–0.38 T; Z −9:9; T = 460 mK | split in 23 interleaved units: atto long (262172), atto express (262173), esbirro (2362), CESGA (11312873); drago never left its queue, its units went to esbirro (2364), CESGA (11313096) and atto (262174) | 274/300 in (atto 262174, units 18–19, running) |

## scan1, preliminary (274 candidates, 2026-11-05)

- Ranking: `scan1_ranking.txt`; lobe 2 alone: `scan1_ranking_L2.txt`; sheets `scan1_top.png`,
  `scan1_top_L2.png`.
- Best overall: g = 15, w = 30 (R_av = 45 nm), μ̃ = 15 meV, τΓ = 3, T_N ≈ 0.23, cost 0.129
  (published 0.174). Lobes 0 and 1 are reproduced (right skew of lobe 1, peak at 0.2 T); g between
  15 and 20 barely matters there.
- **Lobe 2 is not reproduced by the same candidates.** Their I_c is ≈ 0 over 0.25–0.31 T and rises
  towards 0.35 T (right skew); the data peak at 0.27–0.30 T (left skew). Lobe-2 R rms ≈ 0.08 vs
  0.095 for no supercurrent at all.
- The candidates that fit lobe 2 (rms ≈ 0.03: w = 30–35, μ̃ = 10–15 meV, τΓ = 2, T_N ≈ 0.3) lose
  lobe 1: a broad maximum shifted left, the R_0/R_N dip too shallow (0.45 vs 0.25).
- Reading: with the physical Zeeman, no point of this grid reproduces both lobes. The published
  lobe-2 left skew came with physical g ≈ 53.

## Scoring v2 (Carlos, 2026-11-05: weight lobe 2 more, look for the minimal g with a left skew)

cost = rmsR_L1 + 2·rmsR_L2 + 0.5·rel0L. Skew measure `off2`: centroid of the lobe-2 resistance dip,
weight max(1 − R_0/R_N, 0) over 0.245–0.38 T, minus the lobe centre 0.321 T (Φ = 2Φ0 at R_LP = 64 nm);
negative = left skew. **Data: −0.0145 T.** Published curve: −0.0247 T. Per-g tables `scanN_per_g.txt`.

## scan1 final (300 candidates, g 15–20)

Best per g always right-skewed in lobe 2: off2 = +0.0125, +0.0083, +0.0056 T at g = 15, 17.5, 20.
Best overall: g = 20, w = 25, μ̃ = 25, τΓ = 3 (L1 0.078, L2 0.065).

## scan2 (252 candidates: g 20–45 and 17π; w 25–35; μ̃ 10–25; τΓ 2.5–3.5)

Runs: atto long 262175, atto express 262176, esbirro 2365, CESGA 11313960, all saved.
- Fitted radius (w = 25, R_av = 47.5 nm), μ̃ = 25, τΓ = 3, T_N ≈ 0.19: off2 = +0.0056 (g 20), +0.0012 (25),
  −0.0018 (30), −0.0037 (35), −0.0048 (40), −0.0042 (45). The lobe-2 skew turns left at g ≈ 26–27 and
  saturates near −0.005 T; lobe-2 rms 0.065 → 0.044, lobe-1 rms 0.078 → 0.11.
- Check of the sector route at the published parameters (g = 17π, w 25, μ̃ 20, τΓ 2.5): off2 = −0.0254 T
  vs −0.0247 T of the published curve.
- Measured-size left skew (−0.010 to −0.021 T) only in the wide tube (w = 35, R_av = 42.5 nm, μ̃ 10–25,
  τΓ 2.5–3) at any g ≥ 20, with lobe 1 lost (rms 0.17–0.25).
- Between g = 45 (−0.004 T) and 17π (−0.025 T) at the fitted radius the skew jumps; the published
  curve's left skew comes from I_c collapsing above 0.33 T, where the data stay at 5–8 nA.
- Sheet of representatives: `compare.png` (compare.jl).

## gsweep (28 candidates: g 44–56 step 2; w 25; μ̃ 20, 25; τΓ 2.5, 3)

Locates where the strong left skew switches on. atto 262177 (long), 262178 (express), esbirro 2366.
Result (`gsweep_ranking.txt`): the strong left skew needs τΓ = 2.5; at τΓ = 3 off2 stays ≥ −0.005 T up
to g = 56. At τΓ = 2.5, μ̃ = 20: off2 = −0.003 (g 44), −0.006 (46), −0.011 (48), −0.016 (50), −0.022
(52), −0.027 (54), −0.028 (56); μ̃ = 25 within 0.001 T. Smooth onset, no jump: half the measured skew at
g ≈ 46, the measured −0.0145 T at g ≈ 49. Lobe-1 rms 0.12–0.14 throughout (0.078 for the best g = 20,
τΓ = 3 fit). The threshold g drops with τΓ: Zeeman against the coupling-induced CdGM gap (~ZΓ_S).

## alpha (252 candidates: α 0, ±10, ±20, ±30 meV nm; g 15, 17.5, 20; w 25, 30; μ̃ 15–25; τΓ 2.5, 3)

Why: the α = 0 fit needed physical g ≈ 53. In the cylindrical sector Hamiltonian the Rashba term
−(α/r)[m_J − ½σ_z − ½nτ_z + ½(Φ/Φ0)(r²/R_LP²)τ_z]σ_zτ_z has a flux part −(πα R_av/2Φ0)·B·σ_z, a
Zeeman term with g_SO = πα R_av/(μ_BΦ0) ≈ 1.25·α[meV nm] at R_av = 47.5 nm (the orbital+SOC
enhancement of g of Winkler et al., PRL 119, 037701, 2017). g ≈ 17 → 53 would need α ≈ 29 meV nm
with the right sign; measured InAs nanowire α spans 5–32 meV nm. The same term also splits each
sector by ≈ (α/R_av)·m_J, so α ≠ 0 is not just a renormalized g: hence this scan, with both signs
of α since the code's sign of g against InAs's negative g is not settled.
τ calibrated with α included. Smoke test (`smoke_alpha`): α = 0 reproduces the earlier smoke values;
α = 20 changes I_c by ×2.7 at 0.25 T.
Runs: atto 262179 (long, after gsweep), 262180 (express), esbirro 2367, CESGA 11316846.
Outputs: `alpha_ranking.txt`, `alpha_per_alpha.txt` (best per (α, g)), `alpha_top*.png`.
Result (2026-11-05; sheet `compare_alpha.png`, `compare.jl alpha`):
- α < 0 (−10 to −30 meV nm) at g = 15–20 gives a weak left skew, off2 ≈ −0.004 to −0.005 T (a third of
  the data), lobe-2 rms 0.045, lobe-1 rms 0.10–0.12. It saturates already at |α| = 10 and is almost
  independent of g and of |α|: not the Zeeman-like g_SO ≈ 1.25α of the estimate, which would grow with α.
  The m_J-dependent SOC splitting dominates. The lobe-2 current then starts at the lobe edge but stays
  nearly flat (no clear dip).
- α > 0 makes lobe 2 more right-skewed (off2 up to +0.016 T).
- Measured-size left skew at g ≤ 20 only in w = 30, τΓ = 2.5, α = 0 (off2 −0.013 to −0.016 T) with lobe
  1 lost (rms 0.15–0.19).
- Best overall at physical g: α = −10, g = 20, w = 25 (R_av 47.5 nm), μ̃ = 25, τΓ = 3, T_N ≈ 0.175
  (L1 0.097, L2 0.045, off2 −0.004).
- So within the single-radius model SOC does not replace the large g. The measured lobe-2 skew needs
  g_eff ≈ 49 (α = 0, τΓ = 2.5), with lobe-1 rms ≈ 0.13 and I_c collapsing above 0.35 T, where the data
  stay finite. The orbital g enhancement of Winkler et al. comes from SOC between radial subbands, which
  a single-radius model cannot represent: a solid-core (multi-subband) check would be the next test.
