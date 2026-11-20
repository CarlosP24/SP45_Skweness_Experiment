include("plotters/theo_exp_core.jl")

# Teaser image for the PRB Editors' Suggestion (shown as a 200×200 px square thumbnail).
# Top: the full-shell Josephson junction, as a longitudinal cut (SC shell / SM core, as in the
# SP72 sketches) whose two ends open into 3D hexagonal caps (the SP72 / thesis render style)
# that show the CdGM analogs in the core cross-section (Fig. 1b).
# Bottom: the measured dV/dI of Fig. 3c around the 1L and 2L lobes with the theory Ic on top.
# Colour code shared by both panels: green = Ic (supercurrent arrow through the junction and
# the Ic curve), red = B (field arrow right above the map and the map's B axis),
# purple = CdGM analogs.

const color_semi = colorant"#FBF6CB"
const color_super = colorant"#B8C7F9"
const color_core_c = colorant"#D9C54E"        # centre of the core's radial shading
const color_cdgm = colorant"#9B30E0"          # "CdGM analogs" label, as in Fig. 1b
const color_ring1 = colorant"#C23BB5"          # CdGM rings: purple family, kept apart from the green Ic
const color_ring2 = colorant"#5B2BC4"
const color_ic = colorant"#00B140"          # green: contrasts with every part of the seismic map (navy, blue, white, red)
const color_B = :red

hexpts(c, R; sx = 0.82) = [Point2f(c[1] + sx * R * cos(θ), c[2] + R * sin(θ)) for θ in range(0, 2π, 7)[1:6]]

function convexhull(pts)
    P = sort(unique(pts); by = p -> (p[1], p[2]))
    cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
    lower, upper = Point2f[], Point2f[]
    for p in P
        while length(lower) ≥ 2 && cross(lower[end-1], lower[end], p) ≤ 0; pop!(lower); end
        push!(lower, p)
    end
    for p in reverse(P)
        while length(upper) ≥ 2 && cross(upper[end-1], upper[end], p) ≤ 0; pop!(upper); end
        push!(upper, p)
    end
    return vcat(lower[1:end-1], upper[1:end-1])
end

# point inside a convex polygon (vertices in either orientation)
function inconvex(p, poly)
    sg = [(b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1]) for (a, b) in zip(poly, circshift(poly, -1))]
    return all(≥(0), sg) || all(≤(0), sg)
end

# Hexagonal full-shell cap: face at c, prism receding along δ (s = ±1 picks the side)
function hexcap!(ax, c; Ro = 1.35, Ri = 1.0, δ = (0.95, 0.22), s = 1, lw = 2.5)
    d = Point2f(s * δ[1], δ[2])
    face, back = hexpts(c, Ro), hexpts(c, Ro) .+ Ref(d)
    # lateral facet k (vertices k, k+1) faces the viewer if its outward normal points along d
    θv = range(0, 2π, 7)[1:6]
    facet = [cos(θ + π/6) * d[1] + sin(θ + π/6) * d[2] > 0 for θ in θv]
    edge = [facet[k] || facet[mod1(k - 1, 6)] for k in 1:6]     # lateral edge at vertex k

    # translucent prism body and hidden edges
    poly!(ax, convexhull(vcat(face, back)); color = (color_super, 0.55), strokewidth = 0)
    for k in 1:6
        edge[k] || lines!(ax, [face[k], back[k]]; color = (:black, 0.35), linestyle = :dash, linewidth = 1.5)
        facet[k] || lines!(ax, [back[k], back[mod1(k + 1, 6)]]; color = (:black, 0.35), linestyle = :dash, linewidth = 1.5)
    end
    # inner core receding (faint)
    for (k, p) in enumerate(hexpts(c, Ri))
        edge[k] && lines!(ax, [p, p + d]; color = (:black, 0.25), linestyle = :dash, linewidth = 1.2)
    end

    # face: shell, radially shaded core, CdGM rings
    poly!(ax, face; color = color_super, strokewidth = 0)
    inner = hexpts(c, Ri)
    verts = vcat([Point2f(c...)], inner)
    faces = [Makie.GeometryBasics.GLTriangleFace(1, k + 1, mod1(k + 1, 6) + 1) for k in 1:6]
    mesh!(ax, verts, faces; color = vcat([color_core_c], fill(color_semi, 6)), shading = NoShading)

    θs = range(0, 2π, 400)
    ring(r0, a, n) = [Point2f(c[1] + 0.82 * r0 * (1 + a * cos(n * θ)) * cos(θ), c[2] + r0 * (1 + a * cos(n * θ)) * sin(θ)) for θ in θs]
    lines!(ax, ring(0.68Ri, 0.07, 6); color = color_ring1, linewidth = 5)
    lines!(ax, ring(0.38Ri, 0.10, 3); color = color_ring2, linewidth = 5)

    poly!(ax, inner; color = :transparent, strokecolor = :black, strokewidth = lw)
    poly!(ax, face; color = :transparent, strokecolor = :black, strokewidth = lw)

    # visible lateral and back edges
    for k in 1:6
        edge[k] && lines!(ax, [face[k], back[k]]; color = :black, linewidth = lw)
        facet[k] && lines!(ax, [back[k], back[mod1(k + 1, 6)]]; color = :black, linewidth = lw)
    end
    return ring(0.68Ri, 0.07, 6), face
end

function panel_junction(ax; xL = 2.15, xR = 9.85, gap = 0.45, hc = 0.55, ts = 0.32, fs = 34)
    xm = (xL + xR) / 2
    hidedecorations!(ax); hidespines!(ax)

    # longitudinal cut: S1 | bare core | S2
    band!(ax, [xL, xR], -hc, hc; color = color_semi)
    for (a, b) in ((xL, xm - gap), (xm + gap, xR))
        band!(ax, [a, b], hc, hc + ts; color = color_super)
        band!(ax, [a, b], -hc - ts, -hc; color = color_super)
    end

    ringL, faceL = hexcap!(ax, (xL, 0.0); s = -1)
    ringR, faceR = hexcap!(ax, (xR, 0.0); s = 1)

    # supercurrent through the junction (green, like the Ic curve below)
    arrows2d!(ax, [xm - 1.75], [0.0], [3.3], [0.0]; color = color_ic, shaftwidth = 10, tiplength = 30, tipwidth = 34)
    text!(ax, xm + 1.65, 0.0; text = L"I_c", color = color_ic,
        fontsize = fs + 12, align = (:left, :center))

    # B along the wire, right above the map whose horizontal axis it is
    yB = -1.55
    arrows2d!(ax, [xm - 2.6], [yB], [5.0], [0.0]; color = color_B, shaftwidth = 5, tiplength = 22, tipwidth = 22)
    text!(ax, xm + 2.55, yB; text = L"B", color = color_B, fontsize = fs + 4, align = (:left, :center))

    # CdGM analogs label; each arrow leaves a gap after the text and stops short of the cap,
    # aimed at the upper-inner part of the outer ring
    yl = 2.05
    text!(ax, xm, yl; text = "CdGM analogs", color = color_cdgm, fontsize = fs, font = :bold, align = (:center, :center))
    pL = ringL[round(Int, 0.14 * length(ringL))]
    pR = ringR[round(Int, 0.36 * length(ringR))]
    for (p, face, sgn) in ((pL, faceL, -1), (pR, faceR, 1))
        p0 = Point2f(xm + sgn * 2.6, yl - 0.12)
        u = (p - p0) / hypot((p - p0)...)
        t = 0.0
        while !inconvex(p0 + t * u, face); t += 0.01; end   # first point on the cap outline
        p1 = p0 + (t - 0.2) * u
        arrows2d!(ax, [p0[1]], [p0[2]], [p1[1] - p0[1]], [p1[2] - p0[2]]; color = color_cdgm,
            shaftwidth = 4, tiplength = 15, tipwidth = 15)
    end

    xlims!(ax, -0.05, 12.05)
    ylims!(ax, -1.95, 2.45)
    return ax
end

# lobe centroid ∫B Ic dB / ∫Ic dB and maximum over a < B < b
function lobe_stats(B, Ic, a, b)
    idx = findall(x -> a < x < b, B)
    return sum(B[idx] .* Ic[idx]) / sum(Ic[idx]), maximum(Ic[idx])
end

function fig_teaser(exp_data::Exp_Data, name::String; R = 64, Blims = (0.075, 0.395), Ilims = (-46, 46),
        ghost = false, arrows = false, fs = 34, size = (600, 600))
    Trng, _, Is, Brng, Ict, _, _, _, _, _, _, _, dVdImat = process_data(exp_data, "fitted_$(name)")

    Φ0 = 2067.833 #T⋅nm^2
    BΦ(x) = x * Φ0 / (π * R^2)
    Bns = BΦ.([1, 2])
    Bhalfs = BΦ.([0.5, 1.5])

    fig = Figure(; size, fontsize = fs, figure_padding = (6, 10, 4, 6))

    ax = Axis(fig[1, 1]; aspect = DataAspect())
    panel_junction(ax; fs)

    # dV/dI map (Fig. 3c), cropped to the 1L and 2L lobes
    colorrange = (5e2, 1.5e3)
    ax = Axis(fig[2, 1];
        xticks = (Bns, [L"\Phi_0", L"2\Phi_0"]), xticklabelsize = fs, yticksvisible = false, yticklabelsvisible = false,
        spinewidth = 2)
    heatmap!(ax, Trng, Is, dVdImat; colormap = :seismic, colorrange, colorscale = Makie.pseudolog10,
        lowclip = cgrad(:seismic)[1], highclip = cgrad(:seismic)[end], rasterize = 5)

    vlines!(ax, Bns; color = :black, linestyle = :dash, linewidth = 2.5)

    if ghost   # depleted core (no CdGM analogs): symmetric lobes, rescaled to each lobe's height
        Brng_d, Ic_d = ics("depleted_$(name)")
        for (a, b) in ((Bhalfs[1], Bhalfs[2]), (Bhalfs[2], 0.4))
            _, m = lobe_stats(Brng, Ict, a, b)
            idx = findall(x -> a < x < b, Brng_d)
            g = Ic_d[idx] ./ maximum(Ic_d[idx]) .* m
            for sgn in (1, -1)
                lines!(ax, Brng_d[idx], sgn .* g; color = (:white, 0.9), linewidth = 6)
                lines!(ax, Brng_d[idx], sgn .* g; color = :gray15, linewidth = 3.5, linestyle = (:dash, :dense))
            end
        end
    end

    for sgn in (1, -1)
        lines!(ax, Brng, sgn .* Ict; color = color_ic, linewidth = 6)
    end

    # skewness: lobe centre (Φ = nΦ0) → centroid of the Ic lobe
    arrows && for (n, (a, b)) in enumerate(((Bhalfs[1], Bhalfs[2]), (Bhalfs[2], 0.4)))
        Bc, m = lobe_stats(Brng, Ict, a, b)
        y = m + 8
        arrows2d!(ax, [Bns[n]], [y], [Bc - Bns[n]], [0.0]; color = :black, shaftwidth = 6, tiplength = 24, tipwidth = 26)
        scatter!(ax, [Bns[n]], [y]; color = :black, markersize = 13)
    end
    text!(ax, BΦ(2.35), -13; text = L"I_c", color = color_ic,
        fontsize = fs + 12, align = (:center, :top))
    Label(fig[2, 1, Bottom()], L"B"; color = color_B, fontsize = fs + 4, halign = :right, padding = (0, 8, 0, 4))

    xlims!(ax, Blims)
    ylims!(ax, Ilims)

    rowsize!(fig.layout, 1, Relative(0.40))
    rowgap!(fig.layout, 1, 4)
    return fig
end

# fig = fig_teaser(Exp_Data(), "C3")
# save("figures/fig_teaser.pdf", fig)
# save("figures/fig_teaser.png", fig; px_per_unit = 2)   # 1200×1200 px
# fig
