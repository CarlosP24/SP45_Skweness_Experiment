# Poster figures — arXiv:2607.18371

This directory produces the **figures for the A0 conference poster** of
*Probing Synthetic CdGM States Through Critical Current in Full-Shell Nanowire
Josephson Junctions* (Payá et al., 2026).

Scope of this branch: **figures only**. Nothing outside `poster/` changes.
The poster itself is assembled in Affinity Designer from a separate template;
this repo's job is to emit two vector files at exact physical size.

---

## Session start — do this before writing any plotting code

1. Locate the paper's existing figure scripts (search for `Figure`, `CairoMakie`,
   `heatmap`, `fig3`, `LDOS`, `save(`). **Read them.** Do not reimplement from scratch.
2. Locate the data/simulation artefacts they consume (`.jld2`, `.h5`, `.csv`, or a
   Quantica model built in-script). Note whether a figure recomputes or loads.
3. Report what you found and which script maps to which poster slot **before** editing.

The rule for this branch: **re-style, do not re-derive.** The computation that produced
the paper's figures is correct and already fitted. You change typography, geometry,
colour and panel arrangement — not physics.

---

## Non-negotiables

- **Never re-fit.** The parameters are fixed and must match the paper exactly:
  `Rav = 47.5 nm`, `µ = 20 meV`, `Γ_NS = 2.5 Δ₀`, `T_N = 0.24`, `Te = 460 mK`,
  `Δ₀ = 0.2 meV`, `R_LP = 64 nm`, `d = 8 nm`, `ξ_d^L = 55 nm`, `ξ_d^R = 92 nm`,
  `g = 17`, `m* = 0.023 mₑ`, `V_G = 80 V`, `T = 250 mK` for Fig. 2(a) data.
- **Every number that appears on the poster must be traceable to the paper.** If a
  value isn't in the manuscript or its supplement, it does not go on the poster.
- Do not extend field or energy ranges beyond what the paper shows. If a panel looks
  empty at poster aspect, change the layout, not the physics.
- Do not silently drop panels. If something won't fit, say so and propose an option.

---

## The two slots

| Slot | Physical size | CairoMakie `size` (pt) | Content |
|---|---|---|---|
| **Hero** | 771 × 240 mm | `(2185.5, 680.3)` | Paper Fig. 3, rebuilt wide |
| **Supporting** | 373 × 150 mm | `(1057.3, 425.2)` | Paper Fig. 1 (device + model) |

Conversion: `mm × 72/25.4 = pt`. Helper:

```julia
mm2pt(x) = x * 72 / 25.4
const HERO    = (mm2pt(771), mm2pt(240))
const SUPPORT = (mm2pt(373), mm2pt(150))
```

### Why this matters

Export with `pt_per_unit = 1` so **one Makie unit = one PostScript point = the real
printed size**. Then `fontsize = 24` in the figure is genuinely 24 pt on the A0 sheet,
directly comparable to the poster's own type scale (body 32 pt, caption 23 pt). This is
the whole trick — it removes all guesswork about whether labels will be legible.

Never scale a finished figure in Affinity. If the size is wrong, fix `size` here and
re-export.

### Proposed hero layout (confirm before building)

Paper Fig. 3 is a tall 4-row stack; at 3.2:1 a 4-row stack gives 60 mm per panel, far
too flat. Use a **2 × 2 grid**, both columns sharing the `B` axis:

```
┌───────────────────────────┬───────────────────────────┐
│ (a) LDOS, left lead       │ (c) dV/dI zoom + Ic^AH,    │
│                           │     Ic^th overlay          │
├───────────────────────────┼───────────────────────────┤
│ (b) LDOS, right lead      │ (d) R₀/R_N, exp + theory   │
└───────────────────────────┴───────────────────────────┘
        B (T) / Φ-Φ₀ shared          B (T) shared
```

`linkxaxes!` within each column, hide tick labels on the upper row, one x label per
column. Put the `Φ/Φ₀` secondary axis on the top of row 1 only. Label the LP lobes
(0L, 1L, 2L) once, inside panel (a).

### Supporting figure

Paper Fig. 1 is (a) false-colour micrograph over (b) model schematic. At 2.49:1 place
them **side by side**, roughly 180 mm each. The micrograph is a bitmap — check its
native resolution covers ~180 mm at 150 dpi (≈1060 px wide); if not, flag it rather
than upscaling.

---

## Style

Put the theme in `poster/style/poster_theme.jl` and `set_theme!` it at the top of every
figure script. Nothing else sets fonts, sizes or colours.

### Type

Montserrat, matching the poster. Weights: Bold for panel letters, Medium for accents,
Regular for everything else.

```julia
const FONTDIR = joinpath(@__DIR__, "fonts")   # commit the .ttf files here
poster_fonts = (
    regular    = joinpath(FONTDIR, "Montserrat-Regular.ttf"),
    bold       = joinpath(FONTDIR, "Montserrat-Bold.ttf"),
    italic     = joinpath(FONTDIR, "Montserrat-Regular.ttf"),
    bold_italic= joinpath(FONTDIR, "Montserrat-Bold.ttf"),
)
```

| Element | Size (pt) |
|---|---|
| Panel letter (a, b, c, d) | 30, Bold |
| Axis label | 26 |
| Tick label | 22 |
| Colourbar label | 22 |
| In-figure annotation | 26 |
| Legend | 22 |

**20 pt is the absolute floor** — nothing smaller goes in a figure. If a label doesn't
fit at 22 pt, shorten the label.

### Line and mark weights

Data lines 3 pt · axis spines 1.5 pt · ticks 1.5 pt, length 8 · scatter markers 10–14 pt
· grid off. Thin hairlines vanish at two metres; be heavier than feels right on screen.

### Colour

Poster palette — these are the only non-data colours:

```julia
const RED   = colorant"#BF2633"   # highlight / theory curve
const GREEN = colorant"#B5D19D"   # secondary series
const GREY  = colorant"#4A4A4A"   # experiment, axes, all text
const RULE  = colorant"#D8D8D8"
```

Rules:

- **Experiment = `GREY` markers. Theory = `RED` line.** Consistent across every panel,
  so the reader learns it once. This replaces the paper's green/orange/purple/black.
- Red is otherwise reserved for annotation. Do not use it as a data colour anywhere
  else in the same figure — the eye must read red as "look here".
- Heatmaps (LDOS, dV/dI): one perceptually-uniform sequential map per quantity, reused
  wherever that quantity appears. Keep the paper's choice if it is already uniform;
  check it does not put saturated red near the poster's red band. Never a rainbow map.
- Colourbars: label with quantity and unit, ticks at 22 pt, same height as the panel.

### Frames and background

- `backgroundcolor = :transparent` on both `Figure` and `Axis`. The poster supplies the
  white ground.
- **No outer frame.** The template draws a 0.6 mm `#D8D8D8` border around every figure
  slot; a second frame from Makie doubles it.
- No figure title — the poster caption does that job.

---

## Output

```
poster/figures/hero_fig3.pdf          ← goes in the 771 × 240 mm slot
poster/figures/support_fig1.pdf       ← goes in the 373 × 150 mm slot
poster/figures/*.png                  ← 150 dpi previews, for eyeballing only
poster/figures/MANIFEST.md            ← slot → file → source script → data
```

```julia
save("poster/figures/hero_fig3.pdf", fig; pt_per_unit = 1)
save("poster/figures/hero_fig3.png", fig; px_per_unit = 150/72)   # preview
```

Keep heatmaps **vector** (Makie draws them as polygons). Only set `rasterize = 4` on a
heatmap if the PDF exceeds ~30 MB, and say so in `MANIFEST.md` if you do.

Update `MANIFEST.md` on every change — it is how the poster build knows what is current.

---

## Verify before declaring done

Run these, don't assume:

1. **Page size.** `pdfinfo poster/figures/hero_fig3.pdf` → must read `2185.51 x 680.31 pts`
   (= 771 × 240 mm). Any other value means `pt_per_unit` is wrong.
2. **Fonts embedded, no substitution.** `pdffonts poster/figures/hero_fig3.pdf` → every
   row Montserrat, `emb` = yes. A substituted font means the path in the theme is wrong.
3. **Smallest type.** Measure the smallest rendered glyph height; ≥ 20 pt equivalent.
4. **The 1.5 m test.** Crop a representative region at 1:1, print on A4, stand back
   1.5 m. If a tick label is unreadable there, it is unreadable on the poster.
5. **Print-safe.** No pure `#FFFFFF` fills that are meant to be transparent; no stroke
   below 0.3 pt (drops out at press).

---

## Poster context (for judging what a figure has to carry)

The poster's argument, in the order a reader meets it:

> **TL;DR** — A full-shell nanowire in an axial field is a synthetic Abrikosov vortex,
> hosting CdGM analogs. We find their fingerprint in the supercurrent: a critical
> current skewed within every fluxoid lobe.

> **Key takeaway** — The analogs peak at 47.5 nm, the shell at 64 nm. They enclose less
> flux, follow a different Little–Parks period — and that mismatch skews the critical
> current.

So the hero figure must make **two** things visible from three metres, before any text
is read: that Ic is *tilted* within 1L and 2L and *not* in 0L, and that theory tracks it.
If the skew is not obvious at a glance, the figure has failed regardless of how faithful
it is. Consider annotating the skew direction directly on the panel with a red arrow or
a shaded lobe — that annotation is worth more than an extra tick label.

Two honest points that belong on the poster and may need figure support: the 1L
rightward skew is reproducible across devices while the **2L direction is
parameter-dependent** and can go either way; and the abrupt vanishing of Ic at the lobe
edges resembles the destructive LP regime but cannot be it, since destructive regions
would be finite on both sides of the half-integer flux lines.

---

## Conventions

- Branch: `poster/figures`. Commit `.ttf` files and generated PDFs — the poster build
  needs them and they are small.
- One script per slot, named for the slot, not the paper figure number.
- Prefer editing the theme over passing per-call keyword arguments; if two scripts
  disagree about a size or colour, the theme is wrong.
- Ask before changing anything under the paper's own figure directory.
