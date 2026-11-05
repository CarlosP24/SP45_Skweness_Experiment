# Digitized SI data of devices B and C (first estimate)

Stand-in for the lab files of devices B and C, which are not in this repository yet. Recovered from
the Supplemental Material figures of the paper repository (`SI_DEVICE_B.pdf`, `SI_DEVICE_C.pdf`) by
`refit_BC/digitize.jl` with `plots/kitchen/pdf_digitizer.jl` (2026-11-05). Replace them with the
measured files when available; `plots/kitchen/digitized_reader.jl` is the only code that reads them.

| file | source | content |
|---|---|---|
| `devB_dVdI_map.dat` | SI_DEVICE_B panel (a) | dV/dI (Ω), R_s = 5490 Ω subtracted; B = 0:0.004:0.24 T (61) × I = ±750 nA (2069 rows) |
| `devB_R0_line.dat` | SI_DEVICE_B panel (c) | published zero-bias R_0(B), 61 points |
| `devC_dVdI_zoom.dat` | SI_DEVICE_C panel (b) | dV/dI (Ω), R_s = 5350 Ω subtracted; B = 0:0.003:0.48 T (161) × I = −99…101 nA (1769 rows) |
| `devC_dVdI_full.dat` | SI_DEVICE_C panel (a) | same, I = ±800 nA (1769 rows) |
| `devC_R0_line.dat` | SI_DEVICE_C panel (c) | published zero-bias R_0(B), 150 of the 161 points (matplotlib path simplification dropped 11) |
| `digitization_log.txt` | | tick-fit residuals, colour levels, clipping, comparison with the R_0 lines |

Map files: matrix, `#` header; row 1 holds B (T) of each column (first entry NaN), column 1 holds I
(nA) of each row.

## How

- Each map is one embedded RGB image. All its colours are entries of the 256-colour `seismic` table
  of its colourbar (11–17 % of the pixels differ from the colourbar's entry by one unit in one
  channel, a rounding difference; they take the nearest entry). The colour changes at almost every
  pixel: matplotlib interpolated the data before colour-mapping, so the inverted value of a pixel is
  the interpolated measurement there.
- Colour → value: the colourbar image, with its tick marks (vector strokes) fitted linearly. One
  colour level is 7.7–7.9 Ω. Values ≥ 1996 Ω are clipped (device B, 0L switching peaks: 0.3 % of
  the pixels).
- Pixel → (B, I): pixel centres of the image placement, through linear fits of the axis tick marks
  (worst tick residuals 9·10⁻⁶ T and 0.04 nA). Tick values were read off the rendered figure; the
  labels are outlined glyphs.
- Field grid: the R_0 polylines of panel (c) sit on 0.004 T (B) and 0.003 T (C) grids; the map
  columns are interpolated at those fields.

## Checks and caveats

- Device B: the map at I ≈ 0 reproduces the published R_0 line, mean +11 Ω, rms 16 Ω (2 colour
  levels). The map has a current-bias offset of −3.2 nA (zero-bias dip and switching peaks), which
  `refit_BC/target.jl` removes as `plots/fig_gates.jl` does.
- Device C: the zoom map (b) is self-consistent (dip centred at |I| < 1 nA), but the published R_0
  line lies above its I ≈ 0 cut by up to 200 Ω in 2L, rms 76 Ω. The line is closer to the
  full-range map (a) (rms 57 Ω), so it was probably taken from that run, not from (b) as the caption
  says. The fit uses the zoom map, as Fig. 3 uses the low-current run of the main device; the real
  files will settle it.
- The current axis is oversampled (one row per pixel); the step of the original measurement is not
  recoverable from the figure.
