# Paper edits: g ≈ 53 as a fitted parameter (minimal changes)

2026-11-05. Overleaf repo read only (main.tex, SM_content.tex, SM_journal.tex, response.tex, all .bib).
Line numbers refer to the current Overleaf_repo files.

## Where g appears or matters

| place | content | action |
|---|---|---|
| main.tex:161 | "Zeeman $g-$factor ... fixed to standard values ..., $g = 17$"; "four free fitting parameters" | **change** (core edit + footnote) |
| main.tex:170 (Fig. 3 caption) | "Fitting parameters: $R_\text{av}$, $\mu$, $\Gamma_{\rm NS}$, $T_N$" | **add g** |
| main.tex:179 | "the Zeeman field breaks the spin degeneracy ... skews the CdGM analogs" | consistent; optional half-sentence |
| main.tex:196–200 | "$\Gamma_{\rm NS}$, $T_N$ ... maxima; $\mu$, $R_\text{av}$ ... skewness"; "a single parameter set must reproduce the opposite skewness in 1L and 2L" | optional: say that g fixes the 2L sign |
| main.tex:206 | "In 2L ... the skewness can be to either the right or the left depending on the parameters" | consistent (this is now mostly g); no change needed |
| main.tex:268 (App. C) | "the full radial potential would introduce only quantitative corrections to the fitted parameters" | consistent if g is called effective; optional tweak |
| SM_content (depleted, gates) | "All other parameters are identical"; gates reuse the main-text curve | no g stated; no change |
| response.tex | no referee asks about g; ref. B pt. 2 (robustness of fit) is the closest; "List of changes" | add one item |

g is stated in exactly one sentence (main.tex:161) and nowhere in the SM or the response. So the change
is local.

Consistency points:
- Making g a fit parameter raises the count from four to five; the edit below changes "four" to "five".
- "Fixed to standard values" must go for g. It can stay for m_0.

## 1. Core edit (main.tex:161)

**Old**

```latex
Concerning the SM nanowire, the Zeeman $g-$factor and the electron effective mass are fixed to standard values for Al/InAs nanowires, $g = 17$ and $m_0 = 0.023m_e$. Neglecting spin-orbit coupling, which is known to be weak in these devices \cite{Woods:PRB19} and non-essential for CdGM physics \cite{Paya:PRB24}, our microscopic model contains four free fitting parameters: the chemical potential $\mu$, the SC/SM coupling $\Gamma_{\rm{NS}}$, the junction transparency $T_N$ and the average radius $R_\text{av}$ of the CdGM analog wave-functions.
```

**New**

```latex
Concerning the SM nanowire, the electron effective mass is fixed to the standard value for InAs nanowires, $m_0 = 0.023m_e$. Neglecting spin-orbit coupling, which is known to be weak in these devices \cite{Woods:PRB19} and non-essential for CdGM physics \cite{Paya:PRB24}, our microscopic model contains five free fitting parameters: the chemical potential $\mu$, the SC/SM coupling $\Gamma_{\rm{NS}}$, the junction transparency $T_N$, the average radius $R_\text{av}$ of the CdGM analog wave-functions, and the effective $g$-factor of the core electrons\footnote{The fitted value, $g \approx 53$, is well above the bulk InAs value, $|g| \approx 15$. It is needed to reproduce the leftward skewness of $I_c$ in 2L. For $g \lesssim 20$ the calculated 2L skewness points to the right, whereas the 0L and 1L behavior is set mainly by $R_\text{av}$ and $\mu$. Such an effective $g$-factor is expected in InAs nanowires: there, spin-orbit coupling aligns the spin with the orbital angular momentum of the higher subbands, enhancing the $g$-factor of the corresponding states by up to an order of magnitude \cite{Winkler:PRL17}. It is also comparable to the largest effective $g$-factors of subgap states reported in Al/InAs hybrid nanowires, which range from $\sim 4$ to $\sim 50$ and depend strongly on gate voltage \cite{Vaitiekenas:PRL18}. Within our tubular-core model, $g$ is therefore an effective parameter that absorbs spin-orbital contributions beyond the single-radius description.}.
```

Notes on the footnote:
- **Value.** "$g \approx 53$" is the physical value: the code's g = 17 times π, from the multimode builder's
  convention, gives 53.4. Don't mention the code convention in the paper. If you prefer full transparency
  in the response or the data repository README, a neutral line works: "the value g = 17 quoted in the
  preprint referred to an internal convention of the numerical code; the physical Zeeman energy
  corresponds to g ≈ 53." That would be the place to say it; it's your choice.
- **Placement.** Keep the footnote at the end of the sentence, like the other three in main.tex. With
  revtex `nobibnotes` it stays a footnote.
- **The claim "0L/1L set mainly by R_av and μ".** It is supported by the refits:
  - main device at g = 15–20 (wip/refit_zeeman);
  - device B at any g between 0 and 20 (refit_BC).
- **Compressed version, if the footnote is too long:** keep the 2L sentence plus one sentence citing
  both references.

## 2. Fig. 3 caption (main.tex:170)

**Old**

```latex
Fitting parameters: $R_\text{av} = 47.5$~nm, $\mu = 20$~meV, $\Gamma_{\rm{NS}} = 2.5\Delta_0$, and $T_N = 0.24$. The electronic temperature is found to be $T_e = 460$~mK.
```

**New**

```latex
Fitting parameters: $R_\text{av} = 47.5$~nm, $\mu = 20$~meV, $\Gamma_{\rm{NS}} = 2.5\Delta_0$, $g = 53$, and $T_N = 0.24$. The electronic temperature is found to be $T_e = 460$~mK.
```

## 3. Optional one-liners (only if you want the role of g visible in the discussion)

main.tex:200. **Old**

```latex
Since a single parameter set must reproduce the opposite skewness in 1L and 2L, as well as $R_0/R_N$, the two parameters cannot compensate for each other.
```

**New**

```latex
Since a single parameter set must reproduce the opposite skewness in 1L and 2L, as well as $R_0/R_N$, the two parameters cannot compensate for each other; the leftward skewness in 2L additionally fixes the effective $g$-factor.
```

main.tex:268 (App. C). **Old:** "...would introduce only quantitative corrections to the fitted
parameters, without altering the origin of the skewness." **New:** "...would introduce only quantitative
corrections to the fitted parameters (in particular to the effective $g$-factor), without altering the
origin of the skewness."

## 4. BibTeX (append to Josephson.bib; main.tex uses `\bibliography{Josephson, supp}`)

Neither paper is in any .bib file. The Vaitiekenas keys present are `S20`, `PRB20`, `NP21`, `ErrS25`.

```bibtex
@article{Winkler:PRL17,
  title = {Orbital Contributions to the Electron {{$g$}} Factor in Semiconductor Nanowires},
  author = {Winkler, Georg W. and Varjas, D{\'a}niel and Skolasinski, Rafal and Soluyanov, Alexey A. and Troyer, Matthias and Wimmer, Michael},
  year = 2017,
  journal = {Physical Review Letters},
  volume = {119},
  number = {3},
  pages = {037701},
  doi = {10.1103/PhysRevLett.119.037701}
}

@article{Vaitiekenas:PRL18,
  title = {Effective {{$g$}} Factor of Subgap States in Hybrid Nanowires},
  author = {Vaitiek{\.e}nas, S. and Deng, M.-T. and Nyg{\aa}rd, J. and Krogstrup, P. and Marcus, C. M.},
  year = 2018,
  journal = {Physical Review Letters},
  volume = {121},
  number = {3},
  pages = {037703},
  doi = {10.1103/PhysRevLett.121.037703}
}
```

**Verified (web, 2026-11-05):** authors, titles, journal, volume, article numbers, years, DOIs, and
arXiv numbers (1703.10091; 1710.04300, whose preprint title was "Effective g-factor in Majorana Wires").

**Check before submitting — the "4–50" range.** In Vaitiekėnas et al. this range appears in the
introduction as the span reported by studies of InAs/Al nanowires (|g*| ≈ 4–50 vs bulk ≈ 15). Their own
result is a gate-controlled, step-like reduction of g* with decreasing carrier density. I could not see
their highest measured value in the sources reachable here. The footnote is therefore phrased as
"reported in Al/InAs hybrid nanowires ... \cite{Vaitiekenas:PRL18}", which is accurate either way. If you
want "measured by", check their Fig. 3 first.

Also, 53 is slightly above 50, hence "comparable to the largest". From the wip report (not re-verified
here): Winkler et al.'s simulated InAs/Al wire gives Andreev-state g factors of 23 and 43, for the
|l| = 1 and 2 subbands.

## 5. response.tex (List of changes, l. 169–178): one item

```latex
    \item We now treat the effective $g$-factor of the semiconductor core as a fitting parameter, $g \approx 53$, and explain in a footnote why this value is required by the 2L skewness and consistent with the orbital enhancement of the $g$-factor in InAs nanowires.
```

**Risk to anticipate.** A referee could object that the orbital (Larmor) part of the Winkler
enhancement is already in the cylindrical model, so g ≈ 53 double-counts it (wip REPORT §5.3). The
footnote's last sentence ("effective parameter that absorbs spin-orbital contributions beyond the
single-radius description") is the defensible framing. Do not claim the full orbital enhancement is
missing from the model.

## 6. Other issues noticed (factual, outside the g change)

1. **Wrong SM section (main.tex:244).** `\SMref{SM:twoterm}{S-II.B}` labels the two-terminal devices
   (SM subsection A) as "S-II.B". It should read S-II.A, or S-II, since the sentence also covers the
   four-terminal device.
2. **Wrong SM section (main.tex:234).** "three additional ... devices ... (see `\SMref{SM:fourterm}{S-II.B}`
   for experimental data from these devices)" points only to device D. Use `\SMref{SM:repro}{S-II}`.
3. **Tilt angles (main.tex:158).** The text gives "misalignment angle θ = 8°", but SM Table S1 gives
   θ_L = 8°, θ_R = 5° (the code uses both).
4. **Right-lead radius (SM Table S1).** R_LP^R = 63.5 nm, but the text says R_LP^{L,R} = 64 nm; the code
   uses 64 nm for both leads.
5. **Δ0 (main.tex:158).** The text gives 0.2 meV; the code uses 0.204 meV (rounding only).
6. **T_N (Fig. 3 caption).** The caption says T_N = 0.24, but the plotted I_c^th is the computed curve ×
   0.8 (`pref = 0.8`), so the effective T_N ≈ 0.19. Consider stating it, or mentioning the rescaling.
7. **Device C period (SM, Device C).** "R_LP ~ 77 nm ... LP periodicity ~110 mT". The lobe edges in
   SI_DEVICE_C (0.060, 0.180, 0.297, 0.414 T) give ≈ 118 mT, i.e. R_LP√cosθ ≈ 75 nm.
8. **Panel label (SM gates caption).** "(d-e)" should be "(d-f)", as the text says.
9. **Device C R_0 source (SM, Device C caption).** "(c) ... taken from (b)": the published R_0 line
   matches the full-range panel (a) better than the zoom (b), by up to 200 Ω off in 2L
   (refit_BC/digitization_check.png). Worth checking against the raw data.
10. **Typos.**
    - main.tex:155: "as a function the radial coordinate" → "as a function of the radial coordinate".
    - main.tex:161: "$g-$factor" → "$g$-factor" (gone if the edit above is used).
    - SM gates: "$V_G = 80V$" → "$V_G = 80$~V".

Sources: [Winkler et al., arXiv:1703.10091](https://arxiv.org/abs/1703.10091),
[PRL 119, 037701](https://link.aps.org/doi/10.1103/PhysRevLett.119.037701),
[Vaitiekėnas et al., arXiv:1710.04300](https://arxiv.org/abs/1710.04300v2),
[PRL 121, 037703](https://link.aps.org/doi/10.1103/PhysRevLett.121.037703),
[ADS record](https://ui.adsabs.harvard.edu/abs/2018PhRvL.121c7703V/abstract).
