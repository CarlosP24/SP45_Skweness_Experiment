# The Zeeman term in the Paya:26c fit: problem, refit, and the orbital g enhancement

C. Payá (with Claude), 2026-11-05. Concerns the theory curves of *Probing Synthetic CdGM States
Through Critical Current in Full-Shell Nanowire Josephson Junctions* (arXiv:2607.18371, under
review). Working files, scripts and run log: `wip/refit_zeeman/` (`NOTES.md`). Nothing in `src/`
or in the published data has been changed.

## 1. Summary

- The multimode builder `build_cyl_mm` of FullShell.jl, used for every wire in the preprint, puts a
  Zeeman energy π times too large on the diagonal. The published fit, quoted with g = 17, has a
  physical g = 17π ≈ 53.
- Refit with the physical Zeeman (832 candidates; g, w, μ̃, τΓ, and Rashba α): lobes 0 and 1 are
  reproduced at a standard g = 15–20, as well as or better than in the preprint. The left skew of
  the critical current in the second lobe (2L) is not. At α = 0 it needs g ≈ 49; Rashba SOC at
  g = 15–20 gives at most a third of it.
- The orbital enhancement of the g factor in semiconductor nanowires (Winkler *et al.*, PRL 119,
  037701, 2017) does not justify g ≈ 50 in this model. Its orbital part is the Larmor coupling of
  the subband angular momentum, which the cylindrical model already contains exactly. What the
  model may miss is the interband part of the spin-orbit splitting between ±l subbands, and the
  radial structure behind it. Both are captured only by a radially resolved (solid-core) model. That
  is the test that remains.

## 2. The problem

FullShell.jl has two Hamiltonian builders:

| builder | Zeeman term in the code | value |
|---|---|---|
| `build_cyl` (single angular sector `Z`, flux Φ in units of Φ0) | `0.5 g μBΦ0 Φ cosθ / area_LP`, `area_LP = π R_LP²` | ½ g μ_B B cos θ |
| `build_cyl_mm` (all sectors, field B in T) | `0.5 g μBΦ0 πoΦ0 B cosθ`, `πoΦ0 = π/Φ0` | π · ½ g μ_B B cos θ |

with `μBΦ0` = μ_B Φ0 = 119.69 meV nm². In the single-sector builder Φ/Φ0 = πR_LP²B/Φ0 and the
two factors of π cancel. In the multimode builder, recovering μ_B from μ_BΦ0 needs 1/Φ0, but the
code multiplies by π/Φ0. The discrepancy is not the field-versus-flux parametrization: it is an
extra π. Numerically, ½ μ_BΦ0 πoΦ0 = 0.0909 meV/T per unit g, against ½ μ_B = 0.0289 meV/T. The
factor is present in every installed FullShell version, including 1.4.0. It was first noted in
`data/theory/PROVENANCE.md` (2026-09-02), where the m_J-resolved CPR run compensates it with
g = 17π.

**Affected:** every SP45 wire (`Params_mm`), so the LDOS of both leads, I_c(B) of `fitted_C3`, the
depleted-core and gate-dependence referee figures, and the earlier fits. In SP43, only the trivial
fluxoid valve of PRB25a (g = 1 → 3.1). **Unaffected:** PRB25 (g = 0), the topological valve
(g = 0), and all repos without the multimode builder (SP25, SP38, SP72, the thesis figure repos).

The fix is one factor in `HamiltonianBuilder_MM.jl` (`0.5 g μBΦ0 B / Φ0`). It has not been applied.

## 3. Method of the refit

- **Sectors instead of the multimode model.** Both leads have R = 60 nm and d = 8 nm, so they carry
  the same fluxoid number. The junction is cylindrically symmetric, so J_z is conserved through the
  device. The multimode Hamiltonian is then block diagonal in m_J, and Σ_Z of single-sector
  `build_cyl` currents is the same current, now with the physical Zeeman.
  **Validation:** at the published convention (g_1D = 17π), the sector sum reproduces the published
  I_c(B) within 0.5 % of its peak at Φ/Φ0 = 0.25–2. At the published parameters, the lobe-2 skew
  measure (below) is −0.0254 T against −0.0247 T for the published curve.
- **Target:** R_0/R_N(B) and I_c^AH(B) of run `250mK_DEV13_275_16`, the run of Fig. 3, processed as
  published (R_s = 5490 Ω, T_e = 460 mK), plus the dV/dI peak current in lobe 0.

![**Fit target.** Top: dV/dI peak current I_sw (grey) and I_c^AH from R_0/R_N via the
Ambegaokar–Halperin relation at T_e = 460 mK (green), run 275_16. Bottom: R_0/R_N. Lobe 1 is
skewed right (maximum at 0.2 T); lobe 2 has its dip at 0.27–0.30 T, left of its centre at
0.321 T.](target.png){width=75%}

- **Cost:** for each candidate, the overall scale s of I_c is fitted (a transparency rescaling,
  T_N ≈ 0.24 s). The cost is rms[R_th − R_0/R_N] over lobe 1, plus 2 × the same over lobe 2, plus
  0.5 × the relative rms of I_c in lobe 0. R_th is the published Ambegaokar–Halperin conversion.
- **Skew measure off2:** centroid of the lobe-2 resistance dip (weight 1 − R_0/R_N, 0.245–0.38 T)
  minus the lobe centre at Φ = 2Φ0 (0.321 T). Negative means a left skew. **Data: −0.0145 T.**
- Sectors Z = −9…9 (the edge sectors carry < 2·10⁻⁴ of I_c); 76 fields over 0.005–0.38 T;
  T = 460 mK. Runs on atto, esbirro and CESGA.

## 4. Results

**Physical g = 15–20, α = 0** (300 candidates). Lobes 0 and 1 are reproduced, including the right
skew of 1L and the I_c^AH maximum at 0.2 T. Best: g = 20, w = 25 nm (R_av = 47.5 nm), μ̃ = 25 meV,
τΓ = 3, T_N ≈ 0.18, with lobe-1 rms 0.078 against 0.127 for the published curve. Lobe 2 is skewed
**right** in every candidate (off2 = +0.006 to +0.013 T): I_c stays near zero up to about 0.28 T and
rises towards 0.35 T.

**g above 20, α = 0** (scan2 and g sweep, 280 candidates, fitted radius w = 25 nm):

Lobe-2 skew measure off2 (T; data −0.0145) against the physical g:

| g | 20 | 25 | 30 | 35 | 40 | 44 | 46 | 48 | 50 | 52 | 17π | 56 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| τΓ = 3 | +0.006 | +0.001 | −0.002 | −0.004 | −0.005 | −0.005 | −0.004 | −0.002 | −0.001 | −0.001 | – | −0.003 |
| τΓ = 2.5 | | | | | | −0.003 | −0.006 | −0.011 | −0.016 | −0.022 | −0.025 | −0.028 |

(μ̃ = 25 meV in the τΓ = 3 row and 20 meV in the τΓ = 2.5 row; μ̃ = 25 at τΓ = 2.5 is within 0.001 T.)

- At τΓ = 3 the 2L skew turns left at g ≈ 26 and saturates at about −0.005 T.
- At τΓ = 2.5 it grows smoothly from g ≈ 44, reaching the measured value at g ≈ 49. Lobe 1 is then
  worse (rms ≈ 0.13), and I_c collapses above 0.35 T, where the data stay at 5–8 nA. The published
  fit sits in this regime.
- The threshold g drops with τΓ, so the controlling quantity is the Zeeman energy against the
  coupling-induced CdGM gap (∼ZΓ_S), not g alone.

![**Representative candidates at α = 0.** In each panel, top: theory I_c (orange) against I_c^AH
(green) and I_sw (grey); bottom: R_0/R_N, theory (purple) against data (black). Dashed: lobe edges
and the lobe-2 centre. The titles give the parameters, the fitted T_N, the lobe-1 and lobe-2 rms
and off2. From g = 20 (best lobe 1, lobe 2 skewed right) through g = 30–45 at the fitted radius
(weak left skew) to the published parameters at g = 17π (measured-size left skew, worse lobe 1,
I_c collapsing above 0.35 T); bottom right, a wide tube with the left skew but without lobe 1.](compare.png){width=100%}

**Rashba α at g = 15–20** (252 candidates, α = 0, ±10, ±20, ±30 meV nm):

- α < 0 gives a weak left skew, off2 ≈ −0.004 to −0.005 T, with lobe-1 rms 0.10–0.12. It saturates
  already at |α| = 10 meV nm and is almost independent of g and |α|.
- α > 0 makes lobe 2 more right-skewed.
- Best overall at a physical g: α = −10 meV nm, g = 20, w = 25 nm, μ̃ = 25 meV, τΓ = 3, T_N ≈ 0.175.
  Here 2L has current from the lobe edge, but its dip is nearly flat.
- A measured-size left skew at g ≤ 20 appears only for wide tubes (w = 30–35 nm, R_av = 42.5–45 nm,
  τΓ = 2.5), which lose lobe 1 (rms 0.15–0.25).

![**Rashba α at a physical g, against the large-g fits.** Same layout as the previous figure. Top
row: α = 0 at g = 20, then α = −10 (best overall at a physical g) and α = −30: lobe-2 current from
the lobe edge, dip nearly flat, a third of the measured skew. Bottom row: α = +20 (more right
skew), then α = 0 at g = 48 and 50 with τΓ = 2.5, which reach the measured skew.](compare_alpha.png){width=100%}

Full rankings: `scan1_ranking.txt`, `scan2_ranking.txt`, `gsweep_ranking.txt`,
`alpha_ranking.txt`, `alpha_per_alpha.txt`.

## 5. The orbital enhancement of the g factor (Winkler *et al.* 2017)

### 5.1 What the paper shows

G. W. Winkler, D. Varjas, R. Skolasinski, A. A. Soluyanov, M. Troyer and M. Wimmer, *Orbital
contributions to the electron g factor in semiconductor nanowires*, Phys. Rev. Lett. **119**, 037701
(2017), arXiv:1703.10091. They explain why g factors several times the bulk value are measured in
InAs and InSb nanowires, although confinement should reduce g:

- In an approximately cylindrical wire, higher subbands carry orbital angular momentum l. With a
  small effective mass, that gives a large orbital moment, M_o = −(m₀/m\*)μ_B l ê_z (their Eq. 1).
- Strong spin-orbit coupling aligns spin with orbital angular momentum ("L·S-type spin alignment").
  The subbands then form Kramers pairs (l, −l) with opposite spin, split at B = 0 by ΔE_l.
- In an axial field each pair splits with an effective g_{|l|±1/2} = g\* ± 2(m₀/m\*)|l| (their Eq. 2,
  cylindrical, k_z = 0). For InAs, m₀/m\* ≈ 43, so the orbital part alone is ±87|l|.
- The splitting ΔE_l comes from p-like (hole) admixture to the conduction band due to confinement.
  It is computed in an eight-band k·p model and checked with tight binding. It is about 1 meV for
  InAs wires of 40–160 nm.
- The enhancement holds in an axial field below a critical field. Above it, the same-l states bunch
  together and the slope reverts to the bare g\*. For a perpendicular field it almost vanishes (the
  proposed experimental test).
- It appears only at high chemical potential (subbands with |l| ≥ 1 occupied). It is reduced at
  small diameters, optimal between 10 and 100 nm, and survives hexagonal cross-sections and gate
  fields.
- In a simulated InAs wire coupled to Al (40 nm, μ tuned to the |l| = 1 and |l| = 2 subbands), the
  g factor extracted from the slope of the Andreev state is 23 and 43. The bulk InAs value is −14.9.

Experiments on InAs/Al subgap states report effective g factors from about 4 to 50, strongly gate
dependent (S. Vaitiekėnas, M.-T. Deng, J. Nygård, P. Krogstrup and C. M. Marcus, PRL **121**, 037703,
2018). These are slopes of levels, so they include the orbital contribution.

### 5.2 Why it is relevant to this junction

The device sits where Winkler *et al.* predict the largest effect:

- InAs (m\* = 0.023 m₀), in an axial field;
- a radius of 60 nm, inside their 10–100 nm window;
- a chemical potential that fills many angular sectors (μ̃ = 20–25 meV, |m_L| up to 5–6 at
  R_av = 47.5 nm);
- subgap states whose field dependence is what we fit.

A g ≈ 50, read as the slope of a subgap level, is within what this physics produces and what
experiments see.

### 5.3 Why it does not, by itself, justify g ≈ 50 in our model

The orbital part of the enhancement is already in the model. In the cylindrical sector
Hamiltonian, the angular kinetic term (m_l + ½(Φ/Φ0)r²/R_LP² τ_z)²/2m\*r² contains the cross term
m_l·eB/2m\* (ħ = 1). That is exactly −M_o·B of Winkler's Eq. 1: the Larmor coupling of the subband
angular momentum, radius independent. In a full-shell wire it combines with the fluxoid winding
into the reduced flux that disperses and skews the CdGMs (San-Jose *et al.* PRB 107, 155423; thesis
ch3). The fitted g is therefore the spin part alone. Raising it to about 50 to mimic the orbital
contribution would count that contribution twice. The experimental effective g factors of 4–50 are
level slopes, so they compare to our full level dispersion, not to the g parameter.

What remains outside the fitted model is the **spin-orbit splitting ΔE_l** that locks spin to l, and
its field dependence:

- In a conduction-band model of a cylinder, a radial Rashba field gives the term −(α/r)·m_l σ_z τ_z
  in each sector. It splits the ±l Kramers pairs by about 2α|m_l|/R_av, an L·S coupling. This is
  what the α scan tested. α = 10–30 meV nm gives splittings of 0.4–1.3 meV per |m_l| = 1, the scale
  of Winkler's ΔE_l. It produces at most a third of the measured 2L skew. These α are at the upper
  end of what is expected: Woods *et al.* (PRB 99, 161118, 2019) estimate a few meV nm for the
  subbands near the Fermi energy of full-shell InAs/Al wires.
- In Winkler's treatment, ΔE_l comes from the p-like admixture that the confinement induces. It
  depends on the radial shape of each subband's wavefunction and on the electric field at the
  interface. A single-radius model puts every sector on the same circle, R_av, with one constant α.
  It therefore cannot represent sector-dependent SOC splittings, nor the different radii of
  different sectors.

### 5.4 Relation to the solid core

Both missing ingredients are radial. A solid-core junction (thesis ch2–ch4) resolves them within
the conduction-band model:

- **Dome potential and radial Rashba field.** U(r) = U_min + (U_max − U_min)(r/R)², with α(r) =
  −α0 ∂_r U(r). The SOC is strongest at the interface, where the band bending, and with it the
  hole admixture that Winkler identify as the source of ΔE_l, is largest.
- **Sector-dependent radii.** Each sector m_L (and each radial subband m_r) has its own ⟨r²⟩, so its
  own flux mismatch, degeneracy point and L·S splitting. In ch3 this smears the degeneracy points of
  the solid core (`sec:cdgm_smearing`). It changes the skew lobe by lobe, and could produce a 2L left
  skew without a large g.
- **Higher radial subbands.** Winkler's "higher subbands" include radial excitations. A solid core
  can occupy m_r = 1 subbands, which have different orbital content and SOC (thesis ch4 §4.3).
- **Limit.** The interband origin of ΔE_l is beyond a conduction-band model even with a radial
  profile. The solid core captures it only through α0. A full test would use an eight-band k·p core
  (as Winkler *et al.*, or Escribano *et al.*, PRR 2, 033264, 2020).

## 6. Options for the paper

1. **Keep the published curves and state the effective g.** Say that the theory uses a Zeeman
   energy equivalent to g ≈ 50 (17π in the code's convention). Justify it as an effective
   description of spin–orbital effects beyond the single-radius model, citing Winkler *et al.* and
   the measured subgap g factors. Weak point: the orbital part is already in the model, so a referee
   can object that the effective g double counts it.
2. **Refit at a standard g (15–20) and present 2L as qualitative.** Lobes 0 and 1, R_av ≈ 45–47.5 nm
   and the 1L skew are robust. The 2L left skew would then be described as not captured by the
   single-radius model and attributed to radial structure, with a pointer to the solid-core
   calculation. Weak point: Fig. 3(c, d) loses the 2L agreement that the text highlights.
3. **Solid-core junction first, then decide.** Single-sector solid-core junctions (radial sites ×
   19 sectors × 76 fields) at the dome parameters of the thesis, with α0 as a parameter and
   g = 15–20. More expensive than the single-radius scans but feasible on the clusters. If it gives
   the 2L left skew, option 2 becomes a full quantitative fit with a physical g.

In all cases the FullShell bug should be fixed, the affected SP45 data regenerated, and the
referee figures (depleted core, gate dependence) recomputed with the same convention.

## References

- G. W. Winkler *et al.*, Phys. Rev. Lett. **119**, 037701 (2017); arXiv:1703.10091.
- S. Vaitiekėnas *et al.*, Phys. Rev. Lett. **121**, 037703 (2018); arXiv:1710.04300.
- P. San-Jose *et al.*, Phys. Rev. B **107**, 155423 (2023).
- S. D. Escribano, A. Levy Yeyati and E. Prada, Phys. Rev. Research **2**, 033264 (2020).
- B. D. Woods, S. Das Sarma and T. D. Stanescu, Phys. Rev. B **99**, 161118 (2019).
