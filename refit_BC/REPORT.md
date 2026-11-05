# Devices B and C with the physical Zeeman term: m_J-resolved fit of the critical current

C. Payá (with Claude), 2026-11-05. A possible way to handle the Zeeman correction of *Probing
Synthetic CdGM States Through Critical Current in Full-Shell Nanowire Josephson Junctions*
(arXiv:2607.18371): fit the two SI devices instead of relying on the main device's second lobe.
Branch `refit_devices_BC`, folder `refit_BC/` (`NOTES.md` is the run log). The published code and
data are unchanged.

**First estimate.** The lab files of devices B and C are not in the repository. The data used here
were digitized from the SI figures. Both the targets and the fit will move a little with the measured
files.

## Summary

The model is the paper's junction with the Zeeman term corrected: the m_J-independent single-sector
builder summed over m_J, physical g between 0 and 20, no SOC. The fit gives every lobe the same weight
and does not target the lobe-2 skew.

- **Device B is reproduced.**
  - B was measured on 0L and 1L only, on the same chip as the main device.
  - The best fit is g = 7.5, R_av = 43.75 nm, μ̃ = 15 meV, Γ_NS = 3.25 Δ0, T_N = 0.20 and
    T_e = 604 mK.
  - It follows I_c in 0L and I_c and R_0/R_N in 1L: R_0/R_N rms 0.031; 1L skew +17.9 mT against the
    measured +17.3 mT.
  - Any g from 0 to 20 fits nearly as well: in this model the 1L skew is orbital.
- **Device C is only partly reproduced.**
  - C was measured on 0L–3L, on another chip, and its leads were never characterised.
  - 0L and 3L are reproduced.
  - Its flat 1L and its left-skewed 2L are not, at any g ≤ 20. This holds for Γ_NS up to 6 Δ0,
    μ̃ down to 2.5 meV, lead coherence lengths of 30–92 nm, and an alternative R_0 target.
  - Raising g moves the theory's 2L skew away from the data. That is the failure already found for the
    main device's 2L.
- **The 0L–1L fingerprint survives the correction.** The model with a physical g accounts for the
  robust fingerprint of the paper (no skew in 0L, rightward skew in 1L) in device B, as it does in
  the main device (wip). The 2L left skew, seen in A and C, needs physics beyond the single-radius
  model at a physical g.
- **Fig. 3 of the preprint is reproduced** for device B in §4 and for device C in §5.

## 1. Why devices B and C

`wip/refit_zeeman/REPORT.md` found that the multimode builder of FullShell.jl puts a Zeeman energy
π times too large on the diagonal: the published fit (g = 17 in the code) has a physical g ≈ 53. With
the physical Zeeman and g = 15–20, lobes 0 and 1 of the main device (device A) are reproduced, but the
left skew of its second lobe (2L) is not: it needs g ≈ 49, and Rashba SOC gives at most a third of it.

The paper's robust claim is the absence of skew in 0L and the right skew in 1L (referee B, point 3).
Device B was measured only up to the end of 1L. Device C covers lobes 0–4. This fit asks whether the
physical-g model, with SOC switched off, accounts for both devices when the fit does not target the 2L
skew: every lobe n ≥ 1 enters the cost with the same weight and there is no skew term.

## 2. Data

### 2.1 Digitization of the SI figures

`plots/kitchen/pdf_digitizer.jl` reads the vector PDFs of the SI (`SI_DEVICE_B.pdf`, `SI_DEVICE_C.pdf`)
directly: a small interpreter of the page content stream gives every image placement and stroke in
page coordinates.

- **Maps.** Each dV/dI map is one embedded RGB image. All its colours belong to the 256-colour
  `seismic` table of the colourbar, which is a separate embedded image. Inverting that table gives the
  value of each pixel to within one colour level, 7.7–7.9 Ω. The colour changes at almost every pixel,
  so the data were interpolated before colour mapping, and each pixel carries the interpolated
  measurement.
- **Axes.** Fitted to the tick marks, with residuals ≤ 9·10⁻⁶ T and 0.04 nA.
- **Field grid.** Taken from the published R_0(B) lines of panels (c), which are vector polylines:
  0.004 T for B and 0.003 T for C.
- **Clipping.** Values above 2 kΩ are clipped. This affects only device B's 0L switching peaks
  (0.3 % of the pixels).

![**Digitized SI data.** Top: the inverted maps (B panel a, C panels b and a), same colour scale as
the SI. Bottom: published R_0 line (black) against the map at I ≈ 0 (orange). B: rms 16 Ω, two colour
levels. C: the published line sits up to 200 Ω above the zoom map (b) in 2L, and is closer to the
full-range map (a).](figures/digitization_check.png){width=100%}

For device C, the caption says the R_0 line was taken from the zoom map (b), but it matches the
full-range run (a) better. The fit uses the zoom map, just as Fig. 3 uses the low-current run of the
main device (275_16) and not the full-range one of Fig. 2(a). The measured files will settle this.

### 2.2 Processing (as Fig. 3)

`refit_BC/target.jl` runs the published pipeline unchanged on the digitized maps:
`plots/kitchen/data_reader.jl`, `Ic_reader.jl`, and `process_data` of
`plots/plotters/theo_exp_core.jl`. The steps are:

- Savitzky–Golay smoothing along I.
- dV/dI peak current I_sw.
- R_N, the mean dV/dI over I_sw + 20 nA ≤ |I| ≤ I_max.
- R_0/R_N at zero bias.
- I_c^AH = 0.0209 T_e γ_AH(R_0/R_N).

The maps are resampled onto the current step of the main device's run (0.2167 nA), so every
point-count window of the pipeline keeps its width in nA. Two steps had to be added:

- **Current-bias offset.** Device B's map has its zero-bias dip and switching peaks centred at
  −3.2 nA. The offset is removed with `current_offset` of `plots/fig_gates.jl`. Device C's offset is
  −0.05 nA.
- **T_e.** T_e is a fit parameter whose only role is to relate I_c to R_0. As in the paper, it is
  fixed by matching I_c^AH to the resolved dV/dI peak at the 1L maximum. On the main device this
  procedure gives 475 mK, against the published 460 mK. For B it gives 604 mK and for C 292 mK.
  A least-squares match over the whole 1L would give 543 mK on A: away from the maximum, the dV/dI
  peak stops tracking I_c.

![**Fit targets.** Top: dV/dI peak current I_sw (grey) and I_c^AH (green) at the fitted T_e. The
orange band is the 1L window of the T_e fit. Bottom: R_0/R_N. Device A (run 275_16) is shown for
reference.](figures/target.png){width=100%}

Measured skews (centroid of the R_0/R_N dip, minus the lobe centre; + means a right skew):

- **A:** 1L +0.025 T, 2L −0.016 T.
- **B:** 1L +0.017 T.
- **C:** 1L −0.000 T (a broad, flat dip whose minimum is at +0.024 T), 2L −0.011 T, 3L −0.010 T.

Device C also has a left-skewed 2L.

## 3. Model and fit

### 3.1 The m_J-independent builder

The model is the same junction as the paper: tubular core at R_av = R − w/2, Ballistic shell, a0 =
1 nm, semi-infinite leads, T_N calibrated through τ. The difference is that it is built sector by
sector with FullShell's single-sector `build_cyl`, whose Zeeman term is the physical ½gμ_B B cos θ.
Both leads have equal R and d, so J_z is conserved through the junction and Σ_Z of the sector
currents is the multimode current. The code is `src/calculations/Ic_single_mode.jl`; it builds the
wires directly, calibrates τ on the sector-summed normal conductance, and sums the LDOS of the
sectors.

Checks:

- **Device A** (`refit_BC/out/validate.txt`):
  - τ: 0.86242 from the sector sum, against 0.86241 from the multimode model.
  - At the published convention (g_1D = 17π), the sector sum reproduces `data/Js/fitted_C3.jld2`
    within 0.34 % of its peak.
- **Sector cut-offs** (smoke runs): at the largest μ and R_av of each grid, the sectors outside the
  ranges used (|Z| ≤ 9 for B, ≤ 12 for C) carry ≤ 1.8·10⁻⁴ of I_c.

### 3.2 Lead parameters

Devices B and C have no Joule-spectroscopy characterisation. They take device A's d = 8 nm, ξ_d =
55/92 nm, θ = 8°/5° and Δ0 = 0.204 meV, with R_LP set by the measured Little–Parks period:

- **B:** lobe edges at 0.080 and ≈ 0.238 T, giving R_LP = 64.5 nm (R = 60.5 nm). Device B is on the
  same chip as A.
- **C:** edges at 0.060, 0.180, 0.297 and 0.414 T. The period of 0.118 T gives R_LP = 75 nm
  (R = 71 nm). The SM text quotes ≈ 77 nm and ≈ 110 mT; the measured period is 118 mT.

ξ_d and θ of device C are the least certain inputs.

### 3.3 Grid and cost

| | device B | device C |
|---|---|---|
| g (physical) | 0, 5, 10, 15, 20 | 0, 5, 10, 15, 20 |
| R_av (nm) | 52.5, 50, 47.5, 45, 42.5 | 63.5, 61, 58.5, 56, 53.5, 51 |
| μ̃ (meV) | 10, 15, 20, 25, 30 | 10, 15, 20, 25, 30 |
| Γ_NS / Δ0 = τΓ | 2, 2.5, 3, 3.5 | 2, 2.5, 3, 3.5 |
| α | 0 | 0 |
| fields | 0.004–0.24 T, 60 | 0.006–0.42 T, 70 (lobes 0–3) |
| kernel T | 600 mK | 290 mK |
| candidates | 500 | 600 |

For each candidate the cost is

cost = Σ_{n ≥ 1} rms[R_th − R_0/R_N]_{nL} + 0.5 · rms[(s I_c − I_sw)/I_sw(0)]_{0L}

where s is the overall scale of I_c, a transparency rescaling with T_N ≈ 0.24 s, and
R_th = I_0(s I_c/0.0209 T_e/2)^(−2). Every lobe n ≥ 1 has the same weight: the 2L skew is neither
emphasised nor targeted. "T_e fixed" uses the T_e fitted on the data (§2.2), as the paper does; the
kernel temperature of each scan is that value rounded (600 and 290 mK). "T_e free" also fits T_e,
since it only links I_c to R_0. The scans ran on esbirro and atto.

## 4. Device B

### 4.1 Fit

The coarse grid was followed by a refinement around its optimum (81 candidates, half steps;
`B2_fixed_ranking.txt`). The best candidate:

| | g | R_av | μ̃ | Γ_NS | T_N | T_e | rms R_0/R_N (1L) | rel. rms I_sw (0L) | 1L skew, theory / data |
|---|---|---|---|---|---|---|---|---|---|
| **B** | 7.5 | 43.75 nm | 15 meV | 3.25 Δ0 | 0.20 | 604 mK | 0.031 | 0.058 | +17.9 / +17.3 mT |

For comparison, the published main-device fit has R_av = 47.5 nm, μ̃ = 20 meV, Γ_NS = 2.5 Δ0 and
T_N = 0.19 (0.24 × 0.8).

![**Fig. 3 for device B**, same layout as the preprint. (a, b) LDOS at the end of the left and right
semi-infinite leads, summed over m_J sectors with the physical Zeeman. (c) Low-current dV/dI
(digitized), with I_c^AH (green, T_e = 604 mK) and I_c^th (orange). (d) R_0/R_N: data (black) and the
Ambegaokar–Halperin curve of I_c^th (purple). Parameters: g = 7.5, R_av = 43.75 nm, μ̃ = 15 meV,
Γ_NS = 3.25 Δ0, T_N = 0.20, α = 0.](figures/fig_theo_exp_B.png){width=62%}

The theory follows the 0L dV/dI peak down to the lobe edge, the slow rise of I_c^AH across the left
half of 1L, and the right-skewed dip with its minimum at 0.195 T. In this calculation I_c^th vanishes
at the 1L/2L edge.

### 4.2 Role of g

The best candidate at each g of the coarse grid (`B_fixed_per_g.txt`):

| g | 0 | 5 | 10 | 15 | 20 |
|---|---|---|---|---|---|
| cost | 0.085 | 0.074 | 0.067 | 0.078 | 0.075 |
| R_av (nm) / μ̃ (meV) / Γ_NS (Δ0) | 42.5 / 15 / 3.5 | 45 / 15 / 3 | 45 / 15 / 3 | 45 / 15 / 3 | 45 / 20 / 3.5 |
| 1L skew (mT; data +17.3) | +17.5 | +19.2 | +18.9 | +18.1 | +17.7 |

- Device B does not constrain g between 0 and 20. Even g = 0 gives the right 1L skew, because in
  this model the 1L skew is orbital: CdGM analogs at R_av < R_LP see a reduced flux. The Zeeman term
  only refines the step structure.
- The orbital parameters are robust. The top twelve candidates (`figures/B_fixed_top.png`) all have
  R_av = 42.5–47.5 nm, μ̃ = 15–25 meV, Γ_NS = 3–3.5 Δ0 and T_N = 0.17–0.21.
- Letting T_e float lowers the cost only from 0.067 to 0.059, and it does so at T_e ≈ 1.1 K, so T_e
  stays at its data value.

## 5. Device C

### 5.1 Fit with device A's lead parameters

The best candidates on the zoom-map target (`C_fixed_ranking.txt`, `C2_fixed_ranking.txt`):

| scan | g | R_av | μ̃ | Γ_NS | T_N | cost | rms R 1L / 2L / 3L | 0L | skew 1L / 2L / 3L (mT) |
|---|---|---|---|---|---|---|---|---|---|
| data | | | | | | | | | −0.4 / −11.4 / −9.8 |
| coarse | 0 | 56 nm | 15 meV | 3.5 Δ0 | 0.11 | 0.36 | 0.19 / 0.13 / 0.012 | 0.057 | +5.8 / −1.0 / −6.4 |
| extended | 10 | 53.5 nm | 10 meV | 4 Δ0 | 0.12 | 0.36 | 0.15 / 0.12 / 0.044 | 0.072 | +3.3 / −1.7 / −5.5 |

The extended scan (144 candidates) covers τΓ = 4–6 and μ̃ = 2.5–10 meV, past the edges where the
coarse optimum sat.

![**Fig. 3 for device C** with the best candidate overall (extended scan: g = 10, R_av = 53.5 nm,
μ̃ = 10 meV, Γ_NS = 4 Δ0, T_N = 0.12, T_e = 292 mK; lead parameters of device A). Same layout as Fig. 3.
The coarse optimum at g = 0 looks the same.](figures/fig_theo_exp_C.png){width=62%}

C is only partly reproduced, at every g:

- **0L and 3L** are fitted.
- **1L.** The measured I_c^AH is flat across the lobe, 23–26 nA, with R_0/R_N = 0.16–0.22. In every
  candidate, I_c^th starts at 5–10 nA just past the 0L/1L edge and climbs to ≈ 45 nA near 0.14 T.
- **2L.** The left skew (−11 mT) is missed. The theory dip is centred, and the Zeeman term pushes it
  to the right: −1 mT at g = 0, +6 mT at g = 20 (`figures/skew_vs_g.png`). This is the same failure as
  in device A's 2L at a physical g.
- **Targets and T_e do not explain it.** The alternative target (published R_0 line over map (a)'s R_N,
  `C_fixed_Cline_ranking.txt`) has a right-skewed 1L like A and B, but still scores 0.32 with a 1L rms
  of 0.18. Letting T_e float gives 0.34.

![**Lobe skew of the best candidate at each g** (solid) against the data (dashed), with the best cost
(dotted, right axis). Device B's 1L skew is reproduced at every g. Device C's 2L skew moves away from
the data as g grows.](figures/skew_vs_g.png){width=95%}

### 5.2 The lead coherence length

Every C candidate fails the same way: too little supercurrent right after each lobe edge. That points
at the shell gap at the lobe edges, which is set by Λ_edge = ξ_d² Δ0 / (1.76 π R_LP²). Device A's
values at C's R_LP = 75 nm give:

| lead | ξ_d | Λ_edge / Δ0 | shell gap at the 0L/1L edge | same lead on device A (R_LP = 64 nm) |
|---|---|---|---|---|
| left | 55 nm | 0.10 | 0.62 Δ0 | 0.53 Δ0 |
| right | 92 nm | 0.28 | 0.26 Δ0 | 0.09 Δ0 |

On device A, the nearly closed right-lead gap fits the data: I_c^AH is only ≈ 5 nA just past A's 0L/1L
edge. Device C carries 23 nA there, so its leads most likely keep a larger gap at the edges, which
means a shorter ξ_d. To test this, ξ_d of both leads was scanned at g = 10, around the C optimum
(`C3_fixed_ranking.txt`). The scan has 81 candidates: ξ_d = 30, 45, 60 nm (edge gaps 0.84, 0.71,
0.58 Δ0), R_av 53.5–58.5 nm, μ̃ 10–20 meV and Γ_NS 2.5–3.5 Δ0.

**The hypothesis does not hold.** The best ξ_d is 60 nm (cost 0.37), and shorter coherence lengths
make the fit worse:

- **ξ_d = 45 nm:** best cost 0.42.
- **ξ_d = 30 nm:** worse still.
- The 1L rms stays at 0.15–0.18 for every ξ_d.
- The 2L rms grows to 0.19–0.22, because a stronger shell gap deepens and centres the theory's 2L dip.

The missing supercurrent near the lobe edges is therefore not a lead parameter that a shorter ξ_d can
supply. Within the single-radius model, C's flat 1L and its left-skewed 2L stay unexplained at a
physical g, whatever the lead coherence length.

## 6. What this offers for the correction

- **Device B is a clean physical-g replacement for the 0L–1L comparison.**
  - Fig. 3 with device B reproduces I_c in 0L and the right-skewed I_c and R_0/R_N in 1L.
  - It uses g ≤ 20 and no SOC, with orbital parameters close to the published fit.
  - It is the robust fingerprint of the paper: no skew in 0L, rightward skew in 1L (referee B, point 3).
  - B has no 2L data, so the 2L question does not arise.
- **The 2L left skew remains unexplained at a physical g** in both devices that have a measurable
  2L (A and C). The wip analysis traced it to the single-radius model, and a solid-core calculation
  is the open test.
- **Device C is not usable quantitatively yet.**
  - Its leads were never characterised.
  - Its 1L and 2L are reproduced neither with device A's lead parameters nor with a shorter ξ_d
    (§5.2).
  - Its two SI panels disagree on R_0 (§2.1).
  - The lab files and a Joule-spectroscopy fit of its leads come first.
- **Possible structures for the revised Fig. 3:**
  1. Device B as the theory comparison, with device A's Fig. 2 kept as the main measurement. A's
     physical-g refit (wip; 0L and 1L reproduced, 2L qualitative) would go to the SM.
  2. Device A refitted at a physical g in Fig. 3, with device B's fit in the SM as the
     reproducibility check of the 0L–1L agreement.

  Either way, the text quoting "g = 17" has to change: in this model the 1L skew does not need the
  Zeeman term, and any g = 0–20 fits device B.

## 7. Caveats

- **Digitized data.**
  - 7.7–7.9 Ω colour levels and interpolated maps.
  - Device B's 0L switching peaks are clipped at 2 kΩ.
  - Device C's R_0 differs between panels (b) and (c).
- **Assumed lead parameters.** B and C take device A's d, ξ_d, θ and Δ0. This is reasonable for B
  (same chip) and is the main uncertainty for C.
- **T_e.** It comes from the resolved 1L peak (604 mK for B, 292 mK for C), as in the paper, and is
  used as the kernel temperature. B's value is above A's 460 mK although both were measured in the
  same cool-down; T_e here is only the parameter that relates I_c to R_0.
- **Grids.** The grids are coarse in μ̃ (5 meV steps, 2.5 in the refinement). μ̃ changes I_c in
  steps, so nearby values can trade places in the ranking.

## 8. Files and how to reproduce

The branch is `refit_devices_BC`. New code:

- `plots/kitchen/pdf_digitizer.jl`, `plots/kitchen/digitized_reader.jl`
- `src/calculations/Ic_single_mode.jl`
- `refit_BC/*.jl`, `refit_BC/*.sbatch`

Data: `data_exp/digitized/`. No published file was modified. Run from the SP45 root:

```
julia --project=plots refit_BC/digitize.jl            # SI PDFs → data_exp/digitized (needs qpdf)
julia --project=plots refit_BC/target.jl              # experimental targets, T_e
julia --project=.     refit_BC/validate.jl            # single-sector route vs fitted_C3
sbatch -p esbirro -N 8 --ntasks-per-node=32 refit_BC/scan.sbatch B 1/1 es    # grids: B, B2, C, C2, C3
julia --project=plots refit_BC/analyze.jl B2 fixed    # rankings, per-g tables, contact sheets
sbatch refit_BC/final.sbatch B2 fixed 1               # fine I_c and LDOS of the best candidate
julia --project=plots refit_BC/fig_theo_exp_BC.jl fixed 1 B2 C2  # Fig. 3 for B and C
julia --project=plots refit_BC/fig_skew_vs_g.jl
```

The figure script writes `fig_theo_exp_B2.*` and `fig_theo_exp_C2.*`, named after the grids. The
report uses them renamed to `fig_theo_exp_B.*` and `fig_theo_exp_C.*`.

Run log: `refit_BC/NOTES.md`.
