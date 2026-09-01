# Poster theme — the single source of truth for type, colour and geometry.
# `set_theme!(poster_theme)` at the top of every poster figure script; nothing
# else sets fonts, sizes or colours.

using CairoMakie

# ---------------------------------------------------------------- geometry ---
# One Makie unit = one PostScript point, so `fontsize = 24` here is 24 pt on the
# printed A0 sheet. Always export with `pt_per_unit = 1`.
mm2pt(x) = x * 72 / 25.4
const HERO    = (mm2pt(771), mm2pt(240))   # 771 × 240 mm slot
const SUPPORT = (mm2pt(373), mm2pt(150))   # 373 × 150 mm slot
const SUPPORT2 = (0.8SUPPORT[1], 0.75SUPPORT[2])  # 298.4 × 112.5 mm — 80 % × 75 % of it

# ------------------------------------------------------------------- fonts ---
const FONTDIR  = normpath(joinpath(@__DIR__, "..", "fonts"))
const MONT_REG  = joinpath(FONTDIR, "Montserrat-Regular.ttf")
const MONT_MED  = joinpath(FONTDIR, "Montserrat-Medium.ttf")
const MONT_BOLD = joinpath(FONTDIR, "Montserrat-Bold.ttf")

# Montserrat has no Greek and no sub/superscript glyphs. Sub/superscripts are
# solved with rich(), which scales real Montserrat digits — but ω Φ Ω µ have no
# Montserrat glyph at all, in any release. Left to itself Makie would silently
# reach for TeX Gyre Heros and then DejaVu Sans, embedding two extra families,
# so instead one Greek face is committed here and asked for by name. A poster
# PDF therefore holds exactly four fonts: three Montserrat weights, plus
# DejaVu Sans for the handful of Greek glyphs. See MANIFEST.md.
const GREEKFONT = joinpath(FONTDIR, "DejaVuSans.ttf")

poster_fonts = (
    regular     = MONT_REG,
    bold        = MONT_BOLD,
    italic      = MONT_REG,
    bold_italic = MONT_BOLD,
)

# ------------------------------------------------------------------ colours ---
"Wrap a Greek symbol so it comes from the committed Greek face, not a fallback."
gk(s) = rich(s; font = GREEKFONT)

"""
Subscript pinned to the 20 pt floor.

Makie sets sub/superscripts at 0.66 of the parent, which puts a subscript under a
26 pt axis label at 17 pt — below the floor. Pinning to `FS_FLOOR` gives a 0.77
ratio: a shade large for a subscript, and readable from three metres.
"""
sub(s) = subscript(s; fontsize = FS_FLOOR)

hexcolor(s) = parse(Makie.Colors.Colorant, s)

const RED   = hexcolor("#BF2633")   # highlight / theory curve
const GREEN = hexcolor("#B5D19D")   # secondary series
const GREY  = hexcolor("#4A4A4A")   # experiment, axes, all text
const RULE  = hexcolor("#D8D8D8")
const LIGHT = hexcolor("#F2F2F2")   # halos and legend ground on dark heatmaps

# Paper Fig. 3's own data colours, kept at the author's request: inside the hero
# they override the palette's grey-experiment / red-theory rule.
# One pair for both panels, so the reader learns the code once.
const C_EXP   = hexcolor("limegreen")   # experiment: Ic^AH markers and R₀/R_N points
const C_THEO  = hexcolor("orange")      # theory: Ic^th and R₀/R_N curves
const C_HALO  = hexcolor("black")       # the paper's dark under-stroke
const C_ANNOT     = hexcolor("black")    # skew callout, over seismic's pale ground
const C_ANNOT_ALT = hexcolor("#D8D8D8")  # ... and over the dark 1L lobe

# Perceptually uniform, reused wherever the quantity appears.
const CMAP_LDOS = :thermal        # LDOS (a.u.)
const CMAP_DVDI = :davos          # dV/dI — no red anywhere near the poster red

# ------------------------------------------------------------------- sizes ---
const FS_PANEL  = 30              # panel letter, bold
const FS_LABEL  = 26              # axis label
const FS_TICK   = 22              # tick label
const FS_CBAR   = 22              # colourbar label
const FS_ANNOT  = 26              # in-figure annotation
const FS_LEGEND = 22              # legend
const FS_FLOOR  = 20              # absolute floor — nothing renders below this

const LW_DATA  = 3.0
const LW_CURVE = 5.0              # the experiment/theory curves in (c) and (d)
const LW_HALO  = 8.5              # their dark under-stroke
const MS_CURVE = 16               # their markers
const LW_SPINE = 1.5
const LW_RULE  = 2.0              # dashed flux lines
const MS_DATA  = 12

poster_theme = Theme(
    fontsize        = FS_LABEL,
    fonts           = poster_fonts,
    backgroundcolor = :transparent,
    figure_padding  = 12,
    textcolor       = GREY,
    Axis = (
        backgroundcolor  = :transparent,
        xlabelsize = FS_LABEL,  ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK, yticklabelsize = FS_TICK,
        xlabelcolor = GREY, ylabelcolor = GREY,
        xticklabelcolor = GREY, yticklabelcolor = GREY,
        xtickcolor = GREY, ytickcolor = GREY,
        xgridvisible = false, ygridvisible = false,
        spinewidth = LW_SPINE,
        leftspinecolor = GREY, rightspinecolor = GREY,
        topspinecolor = GREY, bottomspinecolor = GREY,
        xtickwidth = LW_SPINE, ytickwidth = LW_SPINE,
        xticksize = 8, yticksize = 8,
        xticklabelpad = 4, yticklabelpad = 4,
    ),
    Lines   = (linewidth = LW_DATA,),
    Scatter = (markersize = MS_DATA,),
    Text    = (fontsize = FS_ANNOT, color = GREY),
    Colorbar = (
        size = 22,
        labelsize = FS_CBAR, ticklabelsize = FS_TICK,
        labelcolor = GREY, ticklabelcolor = GREY,
        spinewidth = LW_SPINE, tickwidth = LW_SPINE, ticksize = 8,
        tickcolor = GREY, spinecolor = GREY,
    ),
    Legend = (
        labelsize = FS_LEGEND, labelcolor = GREY,
        framevisible = false, patchsize = (34, 22),
        padding = (10, 10, 6, 6), rowgap = 2,
    ),
    Label = (color = GREY,),
)

# ---------------------------------------------------------------- helpers ----
"""
Panel letter (a, b, c, d), set outside the panel in its top-left protrusion.

The negative bottom padding is what keeps it outside without cost: a Label in a
protrusion otherwise reserves its own height in the row, and at poster aspect that
is height the panels need.
"""
panel_label!(gp, letter; padding = (0, 0, -30, 0)) =
    Label(gp, letter; font = MONT_BOLD, fontsize = FS_PANEL, color = GREY,
          padding, halign = :left, valign = :bottom)

"Flux quantum in T·nm², and the axial field at n flux quanta through R_LP = 64 nm."
const Φ0_Tnm2 = 2067.833
B_at_flux(n; R = 64) = n * Φ0_Tnm2 / (π * R^2)

"Lobe-centre (integer Φ₀) and lobe-edge (half-integer) fields."
flux_centres(; R = 64) = [B_at_flux(n; R) for n in (1, 2)]
flux_edges(; R = 64)   = [B_at_flux(n; R) for n in (0.5, 1.5)]

"Right-hand edge of the 0L, 1L and 2L lobes, clipped to the plotted field range."
lobe_ends(; R = 64, Bmax = 0.4) = [min(B_at_flux(n; R), Bmax) for n in (0.5, 1.5, 2.5)]

"Dashed Little–Parks guides: solid-dashed at lobe centres, dotted at lobe edges."
function flux_lines!(ax; color = GREY, R = 64)
    vlines!(ax, flux_edges(; R);   color = (color, 0.55), linestyle = :dot,  linewidth = LW_RULE)
    vlines!(ax, flux_centres(; R); color = color,         linestyle = :dash, linewidth = LW_RULE)
    return ax
end

# ------------------------------------------------------------------- output ---
"""
    exact_page!(path, size_pt)

CairoMakie rounds the scene to whole points, so a 771 × 240 mm slot comes out as a
2186 × 680 pt page — 0.17 mm too wide, 0.11 mm too short. Reset the MediaBox to the
exact slot size. qpdf's QDF form is edited in place and repaired with fix-qdf, so every
content stream survives byte for byte; the only thing the width loses is half a point
of empty figure padding.
"""
function exact_page!(path::AbstractString, size_pt::Tuple)
    qpdf, fixqdf = Sys.which("qpdf"), Sys.which("fix-qdf")
    if qpdf === nothing || fixqdf === nothing
        @warn "qpdf/fix-qdf not on PATH — page stays at CairoMakie's rounded size"
        return false
    end
    mktempdir() do dir
        qdf = joinpath(dir, "qdf.pdf")
        run(pipeline(`$qpdf --qdf --object-streams=disable $path $qdf`; stdout = devnull))
        bytes = read(qdf)
        i = first(something(findfirst(codeunits("/MediaBox"), bytes)))
        j = first(something(findnext(codeunits("]"), bytes, i)))
        box = "/MediaBox [ 0 0 $(round(size_pt[1], digits = 4)) $(round(size_pt[2], digits = 4)) ]"
        write(qdf, vcat(bytes[1:i-1], codeunits(box), bytes[j+1:end]))
        fixed = joinpath(dir, "fixed.pdf")
        open(fixed, "w") do io
            run(pipeline(`$fixqdf $qdf`; stdout = io))
        end
        # QDF form is uncompressed — put the streams back the way Cairo had them
        run(pipeline(`$qpdf --stream-data=compress --object-streams=generate $fixed $path`;
                     stdout = devnull))
    end
    return true
end
