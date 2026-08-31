# Supporting slot — 373 × 150 mm. Paper Fig. 1 (device + model), rebuilt side by side.
#
#   ┌──────────────────────┬───────────────────────────────────┐
#   │ (a) false-colour SEM  │ (b) hexagonal cross-section        │
#   │     + contacts, B,    │     + longitudinal junction with   │
#   │     scale, materials  │     U(r), µ and the CdGM analogs   │
#   └──────────────────────┴───────────────────────────────────┘
#
# Paper Fig. 1 is an Affinity Designer document (plots/sketches/fig1.af), not a
# script, so this is a redraw in poster type and colour. The only bitmap is the
# micrograph, extracted from the paper's own fig1.pdf — see MANIFEST.md for its
# resolution limit.
#
# Both panels are drawn in millimetres: one data unit = 1 mm on the printed
# poster, so distances here are literal poster geometry.
#
# Run from anywhere:  julia --project=plots poster/support_fig1.jl

const ROOT = normpath(joinpath(@__DIR__, ".."))
cd(ROOT)

using CairoMakie
using LinearAlgebra: dot

include(joinpath(ROOT, "poster", "style", "poster_theme.jl"))

set_theme!(poster_theme)

# ------------------------------------------------------------------ canvases --
const CANVAS_A = (145.0, 150.0)   # mm
const CANVAS_B = (222.0, 150.0)   # mm
const GAP_AB   = 6.0              # mm  (145 + 222 + 6 = 373)

# ---------------------------------------------------------- micrograph facts --
# Measured off poster/assets/sem_micrograph.jpg: the baked-in scale bar spans
# 83 px for 400 nm, i.e. 4.82 nm/px. Rows 226:285 hold that bar and its Arial
# label, so they are cropped away and the bar is redrawn in Montserrat.
const SEM_PATH    = joinpath(ROOT, "poster", "assets", "sem_micrograph.jpg")
const SEM_KEEP    = 225           # rows kept from the top
const SEM_NM_PX   = 400 / 83      # nm per pixel

# false colours, sampled from the paper's own Fig. 1 legend chips
const C_AL   = hexcolor("#b8c7f9")
const C_INAS = hexcolor("#fbf6cb")
const C_CORE = hexcolor("#d3c653")
const C_TIAL = hexcolor("#e2dcff")

# The one geometric fact the schematic asserts, from the paper's parameters:
# the analogs peak at R_av = 47.5 nm inside an R = 60 nm core, and the shell
# where they are classically allowed reaches in to about 0.68 R.
const RAV    = 47.5 / 60
const RALLOW = 0.68

blank_axis(gp, w, h) = begin
    ax = Axis(gp; aspect = DataAspect(), backgroundcolor = :transparent)
    hidedecorations!(ax); hidespines!(ax)
    limits!(ax, 0, w, 0, h)
    ax
end

"Regular hexagon, flat top and bottom, circumradius R about centre c."
hexagon(c, R) = [Point2f(c[1] + R*cos(k*π/3), c[2] + R*sin(k*π/3)) for k in 0:5]

"Closed curve at radius r0 with the six-fold modulation of the hexagonal core."
function hexring(c, r0; amp = 0.07, n = 240)
    θ = range(0, 2π, n)
    [Point2f(c[1] + r0*(1 + amp*cos(6t))*cos(t), c[2] + r0*(1 + amp*cos(6t))*sin(t)) for t in θ]
end

# ------------------------------------------------------------- (a) micrograph --
# A 2.7:1 image in a 145 × 150 mm column cannot fill the column's height, so the
# whole group — gate above, image, materials below — is centred as one block.
function panel_device!(ax)
    img = rotr90(Makie.FileIO.load(SEM_PATH)[1:SEM_KEEP, :])
    npx, npy = size(img)

    x0, x1 = 12.0, 136.0                       # 124 mm wide
    y0     = 53.5
    y1     = y0 + (x1 - x0) * npy / npx        # keep the pixel aspect
    image!(ax, (x0, x1), (y0, y1), img; interpolate = true)

    mmpx = (x1 - x0) / npx                     # mm per source pixel
    xat(px) = x0 + px * mmpx
    yat(px) = y1 - px * mmpx
    ywire = yat(145)                           # the nanowire's row in the image

    # contacts: current in on the left, voltage probe on the right
    text!(ax, 2, ywire; text = "I", align = (:left, :center), font = MONT_MED)
    arrows2d!(ax, [10.0], [ywire], [x0 + 9], [ywire]; argmode = :endpoint,
              color = GREY, shaftwidth = 2, tipwidth = 10, tiplength = 9)
    text!(ax, 142, ywire; text = "V", align = (:right, :center), font = MONT_MED)
    lines!(ax, [x1 - 9, 138.5], [ywire, ywire]; color = GREY, linewidth = 2)

    # gate, above the middle finger
    text!(ax, xat(310), 105.5; text = rich("V", subscript("G")),
          align = (:center, :bottom), font = MONT_MED)
    lines!(ax, [xat(310), xat(310)], [103.5, yat(45)]; color = GREY, linewidth = 2)

    # ground reference, in the same relative spot as the paper's Fig. 1
    gx, gy = xat(445), 100.0
    lines!(ax, [gx, gx], [gy, gy + 5]; color = GREY, linewidth = 2)
    for (i, hw) in enumerate((5.0, 3.2, 1.6))
        lines!(ax, [gx - hw, gx + hw], [gy - 2.4(i-1), gy - 2.4(i-1)];
               color = GREY, linewidth = 2)
    end

    # axial field and scale, in the dark band below the wire
    yB = y0 + 7
    text!(ax, xat(300), yB; text = "B", align = (:right, :center), color = LIGHT,
          font = MONT_MED)
    arrows2d!(ax, [xat(325)], [yB], [xat(465)], [yB]; argmode = :endpoint,
              color = LIGHT, shaftwidth = 2.5, tipwidth = 11, tiplength = 10)

    bar = 400 / SEM_NM_PX * mmpx               # 400 nm at the measured 4.82 nm/px
    bx  = xat(115)
    lines!(ax, [bx, bx + bar], [y0 + 3.5], [y0 + 3.5]; color = LIGHT, linewidth = 3.5)
    text!(ax, bx + bar/2, y0 + 5; text = "400 nm", align = (:center, :bottom),
          color = LIGHT, fontsize = FS_TICK)

    # materials
    for (i, (col, lbl)) in enumerate(((C_AL, "Al"), (C_INAS, "InAs"), (C_TIAL, "Ti/Al")))
        cx = 14.0 + (i - 1) * 44
        poly!(ax, Rect2f(cx, 36, 7, 7); color = col, strokecolor = GREY, strokewidth = 1.2)
        text!(ax, cx + 10, 39.5; text = lbl, align = (:left, :center), fontsize = FS_TICK)
    end
    return ax
end

# ------------------------------------------------------------------ (b) model --
function panel_model!(ax)
    # ---- hexagonal cross-section, extruded a little so it reads as a wire
    c        = (58.0, 81.0)
    Rsh, Rco = 42.0, 29.0
    off      = Point2f(-13, -11)
    outer, core = hexagon(c, Rsh), hexagon(c, Rco)
    back = [p + off for p in outer]
    poly!(ax, back; color = hexcolor("#93a6d6"), strokecolor = GREY, strokewidth = 1.2)
    for k in 1:6                               # only the faces turned towards us
        p, q = outer[k], outer[mod1(k + 1, 6)]
        n = Point2f(q[2] - p[2], p[1] - q[1])  # outward normal
        dot(n, off) > 0 && poly!(ax, [p, q, q + off, p + off];
                                 color = hexcolor("#93a6d6"), strokecolor = GREY,
                                 strokewidth = 1.2)
    end
    poly!(ax, outer; color = C_AL,   strokecolor = GREY, strokewidth = 1.5)
    poly!(ax, core;  color = C_INAS, strokecolor = GREY, strokewidth = 1.5)

    # green marks the shell where the states are allowed (U < µ); red is where
    # the analogs peak — R_av = 47.5 nm of the R = 60 nm core, i.e. 0.79 R.
    band = vcat(hexring(c, Rco; amp = 0.0), reverse(hexring(c, RALLOW * Rco)))
    poly!(ax, band; color = (GREEN, 0.9), strokecolor = :transparent)
    lines!(ax, hexring(c, RAV * Rco); color = RED, linewidth = LW_DATA)

    # ---- longitudinal junction
    xL, xR   = 112.0, 216.0
    yc       = 81.0
    hcore    = 26.0
    hshell   = 15.0
    xj0, xj1 = 180.0, 196.0                    # bare, gated weak link
    poly!(ax, Rect2f(xL, yc - hcore, xR - xL, 2hcore); color = C_INAS,
          strokecolor = GREY, strokewidth = 1.2)
    for (a, b) in ((xL, xj0), (xj1, xR)), s in (+1, -1)
        y = s > 0 ? yc + hcore : yc - hcore - hshell
        poly!(ax, Rect2f(a, y, b - a, hshell); color = C_AL,
              strokecolor = GREY, strokewidth = 1.2)
    end
    text!(ax, xR - 4, yc - hcore - hshell/2; text = "Al", align = (:right, :center),
          fontsize = FS_ANNOT)
    text!(ax, xR - 4, yc - 14; text = "InAs", align = (:right, :center), fontsize = FS_ANNOT)

    # gate potential over the weak link
    xg = range(xj0 - 14, xj1 + 14, 160)
    lines!(ax, xg, yc .+ hcore .+ 2 .+ 13 .* exp.(-((xg .- (xj0 + xj1)/2) ./ 7.0).^2);
           color = (GREY, 0.8), linewidth = 2.5)
    text!(ax, xj1 + 12, yc + hcore + 11; text = rich("V", subscript("G")),
          align = (:left, :center), fontsize = FS_ANNOT)

    # radial potential: a hill on the axis, so the states pile up near the surface
    r      = range(-hcore, hcore, 300)
    σ      = 20.0
    U      = exp.(-(r ./ σ).^2)
    xU, wU = 116.0, 34.0
    rav    = RAV * hcore                       # R_av / R = 47.5 / 60
    rall   = RALLOW * hcore                    # classical turning point, U = µ
    μ̂      = exp(-(rall / σ)^2)
    xμ     = xU + wU*μ̂

    shell = [Point2f(xμ, yc + rr) for rr in range(rall, hcore, 40)]
    append!(shell, [Point2f(xU + wU*exp(-(rr/σ)^2), yc + rr) for rr in range(hcore, rall, 40)])
    poly!(ax, shell; color = (GREEN, 0.9), strokecolor = :transparent)
    poly!(ax, [Point2f(p[1], 2yc - p[2]) for p in shell]; color = (GREEN, 0.9),
          strokecolor = :transparent)

    lines!(ax, xU .+ wU .* U, yc .+ r; color = GREY, linewidth = LW_DATA)
    lines!(ax, [xμ, xμ], [yc - hcore - 26, yc + hcore]; color = GREY,
           linewidth = 2, linestyle = :dash)
    text!(ax, xμ, yc - hcore - 28; text = "µ", align = (:center, :top), fontsize = FS_ANNOT)
    text!(ax, xU + wU + 3, yc + 12; text = "U(r)", align = (:left, :center),
          fontsize = FS_ANNOT)

    ψ = exp.(-((abs.(r) .- rav) ./ 3.2).^2)
    lines!(ax, xU .+ 30 .* ψ, yc .+ r; color = RED, linewidth = LW_DATA)

    # axial field
    text!(ax, 118, 136; text = "B", align = (:right, :center), font = MONT_MED)
    arrows2d!(ax, [124.0], [136.0], [xR], [136.0]; argmode = :endpoint,
              color = GREY, shaftwidth = 2.5, tipwidth = 11, tiplength = 10)

    # the callout the poster is built around
    text!(ax, 78, 14; text = "CdGM analogs", align = (:center, :bottom),
          color = RED, font = MONT_BOLD)
    arrows2d!(ax, [66.0], [23.0], [c[1] + 0.58Rco], [c[2] - 0.72Rco]; argmode = :endpoint,
              color = RED, shaftwidth = 2, tipwidth = 9, tiplength = 8)
    arrows2d!(ax, [96.0], [23.0], [xU + 22], [yc - rav]; argmode = :endpoint,
              color = RED, shaftwidth = 2, tipwidth = 9, tiplength = 8)
    return ax
end

function support_fig1()
    fig = Figure(size = SUPPORT, figure_padding = 0)
    axA = blank_axis(fig[1, 1], CANVAS_A...)
    axB = blank_axis(fig[1, 2], CANVAS_B...)
    panel_device!(axA)
    panel_model!(axB)

    colsize!(fig.layout, 1, Relative(CANVAS_A[1] / (CANVAS_A[1] + CANVAS_B[1] + GAP_AB)))
    colgap!(fig.layout, 1, GAP_AB * 72 / 25.4)

    panel_letter!(axA, "a")
    panel_letter!(axB, "b")
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    fig = support_fig1()
    mkpath("poster/figures")
    save("poster/figures/support_fig1.pdf", fig; pt_per_unit = 1)
    save("poster/figures/support_fig1.png", fig; px_per_unit = 150/72)
    @info "wrote poster/figures/support_fig1.{pdf,png}"
end
