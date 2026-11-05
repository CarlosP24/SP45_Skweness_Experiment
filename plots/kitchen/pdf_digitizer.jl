# Digitizer for the vector + raster figure PDFs of the Supplemental Material (devices B and C).
#
# Stop-gap until the lab files of devices B and C are available: it recovers, from the published
# SI figures, the same arrays the lab files give for the main device (dV/dI(B, I) after the R_s
# subtraction, and the zero-bias R_0(B) curve), so that the experimental pipeline of
# plots/kitchen/data_reader.jl and Ic_reader.jl can run on them unchanged.
#
# How the figures are built (checked on SI_DEVICE_B.pdf and SI_DEVICE_C.pdf, Inkscape 1.3.2/cairo):
#   * each dV/dI map is one embedded RGB image (≈ 2000 × 2000 px). Every pixel colour is an exact
#     entry of the 256-colour `seismic` lookup table, and the colour changes at almost every pixel:
#     the data were interpolated *before* colour mapping (matplotlib's default interpolation stage).
#     Inverting the lookup table therefore returns the interpolated data value at each pixel, to
#     within one colour level, (vmax − vmin)/256 = 7.8 Ω for the 0–2 kΩ colourbars;
#   * the colourbars are separate embedded images whose rows sweep the same lookup table, with tick
#     marks as vector strokes: they give the colour → value map;
#   * axes, ticks and the R_0(B) line are vector strokes; tick labels are outlined glyphs, so the
#     tick *values* are given by hand (read off the rendered figure) and their *positions* are fitted.
#
# Requires `qpdf` (to write an uncompressed copy without object streams that this parser can read).
# No dependency beyond Base.

module PDFDigitizer

export digitize_pdf, PathRec, ImageRec, image_pixels, fit_axis, ticks_below, ticks_left,
       ticks_right, frame_of, colorbar_lut, invert_image, pixel_centres, image_extent, strokes_in,
       polyline

# ---------------------------------------------------------------------------------------------
# PDF objects (qpdf --qdf output)
# ---------------------------------------------------------------------------------------------

struct Obj
    dict::String
    stream::Union{Nothing,Vector{UInt8}}
end

function load_qdf(path)
    bytes = read(path)
    s = String(copy(bytes))
    objs = Dict{Int,Obj}()
    for m in eachmatch(r"\n(\d+) 0 obj\n", s)
        id = parse(Int, m.captures[1])
        start = m.offset + length(m.match)
        body_end = first(findnext("\nendobj", s, start))
        k = findnext("\nstream\n", s, start)
        if k !== nothing && first(k) < body_end
            objs[id] = Obj(s[start:first(k)-1], bytes[last(k)+1:first(findprev("endstream", s, body_end))-1])
        else
            objs[id] = Obj(s[start:body_end], nothing)
        end
    end
    return objs
end

ref(dict, key) = (m = match(Regex("/" * key * " (\\d+) 0 R"), dict); m === nothing ? nothing : parse(Int, m.captures[1]))

function stream_data(objs, id)
    o = objs[id]
    L = ref(o.dict, "Length")
    n = L === nothing ? parse(Int, match(r"/Length (\d+)", o.dict).captures[1]) : parse(Int, strip(objs[L].dict))
    return o.stream[1:n]
end

function resource_map(objs, dict, sub)
    out = Dict{String,Int}()
    r = ref(dict, "Resources")
    rd = r === nothing ? dict : objs[r].dict
    m = match(Regex("/" * sub * "\\s*<<(.*?)>>", "s"), rd)
    body = if m === nothing
        r2 = ref(rd, sub)
        r2 === nothing && return out
        objs[r2].dict
    else
        m.captures[1]
    end
    for mm in eachmatch(r"/(\w+)\s+(\d+) 0 R", body)
        out[mm.captures[1]] = parse(Int, mm.captures[2])
    end
    return out
end

function numarray(dict, key)
    m = match(Regex("/" * key * "\\s*\\[([^\\]]*)\\]"), dict)
    return m === nothing ? nothing : parse.(Float64, split(strip(m.captures[1])))
end

# ---------------------------------------------------------------------------------------------
# Content-stream interpreter: graphics state, paths and image placements in page coordinates
# (PDF user space of the page: origin bottom-left, y up, 1 unit = 1 pt)
# ---------------------------------------------------------------------------------------------

const Mat = NTuple{6,Float64}
mul(m1::Mat, m2::Mat) = (m1[1]*m2[1] + m1[2]*m2[3], m1[1]*m2[2] + m1[2]*m2[4],
                         m1[3]*m2[1] + m1[4]*m2[3], m1[3]*m2[2] + m1[4]*m2[4],
                         m1[5]*m2[1] + m1[6]*m2[3] + m2[5], m1[5]*m2[2] + m1[6]*m2[4] + m2[6])
apply(m::Mat, x, y) = (m[1]*x + m[3]*y + m[5], m[2]*x + m[4]*y + m[6])

mutable struct GState
    ctm::Mat
    stroke::NTuple{3,Float64}
    fill::NTuple{3,Float64}
    lw::Float64
    dash::Vector{Float64}
end
Base.copy(g::GState) = GState(g.ctm, g.stroke, g.fill, g.lw, copy(g.dash))

"A painted path: `kind` is :stroke or :fill; points in page coordinates; line width and dash in pt."
struct PathRec
    kind::Symbol
    subpaths::Vector{Vector{NTuple{2,Float64}}}
    color::NTuple{3,Float64}
    lw::Float64
    dash::Vector{Float64}
end

"An image draw: object id, size, and the CTM mapping the unit square onto the page."
struct ImageRec
    id::Int
    W::Int
    H::Int
    ctm::Mat
end

function tokenize(s::AbstractString)
    toks = String[]
    i, n = firstindex(s), lastindex(s)
    while i <= n
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
        elseif c == '%'
            j = findnext('\n', s, i); i = j === nothing ? n + 1 : j
        elseif c == '('
            depth, j = 0, i
            while j <= n
                cj = s[j]
                if cj == '\\'
                    j = nextind(s, j)
                elseif cj == '('
                    depth += 1
                elseif cj == ')'
                    depth -= 1
                    depth == 0 && break
                end
                j = nextind(s, j)
            end
            push!(toks, "(str)"); i = nextind(s, j)
        elseif c == '[' || c == ']'
            push!(toks, string(c)); i = nextind(s, i)
        elseif c == '<' && i < n && s[nextind(s, i)] == '<'
            push!(toks, "<<"); i = nextind(s, nextind(s, i))
        elseif c == '>' && i < n && s[nextind(s, i)] == '>'
            push!(toks, ">>"); i = nextind(s, nextind(s, i))
        elseif c == '<'
            j = findnext('>', s, i); push!(toks, "(hex)"); i = nextind(s, j)
        else
            j = i
            while j <= n && !isspace(s[j]) && !(s[j] in "[]()<>%") && !(j > i && s[j] == '/')
                j = nextind(s, j)
            end
            push!(toks, s[i:prevind(s, j)]); i = j
        end
    end
    return toks
end

function interpret!(paths, images, objs, content::String, resdict::String, gs::GState)
    xobj = resource_map(objs, resdict, "XObject")
    stack, gstack = Any[], GState[]
    cur = Vector{Vector{NTuple{2,Float64}}}()
    sub = NTuple{2,Float64}[]
    inarray, arr = false, Float64[]
    scale(g) = sqrt(abs(g.ctm[1]*g.ctm[4] - g.ctm[2]*g.ctm[3]))
    endsub!() = (isempty(sub) || push!(cur, sub); sub = NTuple{2,Float64}[])
    function paint!(kinds...)
        endsub!()
        for k in kinds
            push!(paths, k == :stroke ? PathRec(:stroke, cur, gs.stroke, gs.lw * scale(gs), gs.dash .* scale(gs)) :
                                        PathRec(:fill, cur, gs.fill, 0.0, Float64[]))
        end
        cur = Vector{Vector{NTuple{2,Float64}}}()
    end
    for t in tokenize(content)
        if t == "["
            inarray, arr = true, Float64[]; continue
        elseif t == "]"
            inarray = false; push!(stack, copy(arr)); continue
        elseif inarray
            (x = tryparse(Float64, t)) === nothing || push!(arr, x); continue
        end
        x = tryparse(Float64, t)
        if x !== nothing || startswith(t, "/") || t in ("(str)", "(hex)", "<<", ">>")
            push!(stack, x === nothing ? t : x); continue
        end
        v = Float64[y for y in stack if y isa Float64]
        if t == "q"
            push!(gstack, copy(gs))
        elseif t == "Q"
            gs = pop!(gstack)
        elseif t == "cm"
            gs.ctm = mul(Mat(v[end-5:end]), gs.ctm)
        elseif t == "RG"
            gs.stroke = Tuple(v[end-2:end])
        elseif t == "rg"
            gs.fill = Tuple(v[end-2:end])
        elseif t == "G"
            gs.stroke = (v[end], v[end], v[end])
        elseif t == "g"
            gs.fill = (v[end], v[end], v[end])
        elseif t == "w"
            gs.lw = v[end]
        elseif t == "d"
            a = [y for y in stack if y isa Vector{Float64}]
            gs.dash = isempty(a) ? Float64[] : last(a)
        elseif t == "m"
            endsub!(); push!(sub, apply(gs.ctm, v[end-1], v[end]))
        elseif t == "l"
            push!(sub, apply(gs.ctm, v[end-1], v[end]))
        elseif t == "c"
            push!(sub, apply(gs.ctm, v[end-1], v[end]))      # end point only: curves are glyphs
        elseif t == "v" || t == "y"
            push!(sub, apply(gs.ctm, v[end-1], v[end]))
        elseif t == "h"
            isempty(sub) || push!(sub, sub[1])
        elseif t == "re"
            endsub!()
            x0, y0, w, h = v[end-3:end]
            push!(cur, [apply(gs.ctm, x0, y0), apply(gs.ctm, x0 + w, y0), apply(gs.ctm, x0 + w, y0 + h),
                        apply(gs.ctm, x0, y0 + h), apply(gs.ctm, x0, y0)])
        elseif t in ("S", "s")
            paint!(:stroke)
        elseif t in ("f", "F", "f*")
            paint!(:fill)
        elseif t in ("B", "B*", "b", "b*")
            paint!(:fill, :stroke)
        elseif t == "n"
            endsub!(); cur = Vector{Vector{NTuple{2,Float64}}}()
        elseif t == "Do"
            name = String(last([y for y in stack if y isa String && startswith(y, "/")]))[2:end]
            id = xobj[name]
            o = objs[id]
            if occursin("/Subtype /Image", o.dict)
                W = parse(Int, match(r"/Width (\d+)", o.dict).captures[1])
                H = parse(Int, match(r"/Height (\d+)", o.dict).captures[1])
                push!(images, ImageRec(id, W, H, gs.ctm))
            elseif occursin("/Subtype /Form", o.dict)
                M = numarray(o.dict, "Matrix")
                g2 = copy(gs)
                M === nothing || (g2.ctm = mul(Mat(M), g2.ctm))
                interpret!(paths, images, objs, String(stream_data(objs, id)), o.dict, g2)
            end
        end
        empty!(stack)
    end
    return paths, images
end

"""
    digitize_pdf(pdf) -> (objs, paths, images)

Normalise `pdf` with qpdf, then interpret its (single) page: every painted path and every image
placement, in page coordinates.
"""
function digitize_pdf(pdf::AbstractString)
    qdf = tempname() * ".pdf"
    run(`qpdf --qdf --object-streams=disable $pdf $qdf`)
    objs = load_qdf(qdf)
    rm(qdf; force = true)
    pid = first(k for (k, o) in objs if occursin(r"/Type /Page\b", o.dict))
    pdict = objs[pid].dict
    content = String(stream_data(objs, ref(pdict, "Contents")))
    paths, images = PathRec[], ImageRec[]
    gs = GState((1.0, 0.0, 0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 1.0, Float64[])
    interpret!(paths, images, objs, content, pdict, gs)
    return objs, paths, images
end

"Pixels of an image object as an H × W × 3 UInt8 array; row 1 is the top of the image."
function image_pixels(objs, id)
    o = objs[id]
    W = parse(Int, match(r"/Width (\d+)", o.dict).captures[1])
    H = parse(Int, match(r"/Height (\d+)", o.dict).captures[1])
    occursin("/Filter", o.dict) && error("image $id is still encoded: $(o.dict)")
    occursin(r"/BitsPerComponent 8", o.dict) && occursin(r"/ColorSpace /DeviceRGB", o.dict) ||
        error("image $id is not 8-bit RGB")
    d = stream_data(objs, id)
    length(d) == 3W * H || error("image $id: $(length(d)) bytes for $W × $H RGB")
    return permutedims(reshape(d, 3, W, H), (3, 2, 1))
end

"Page-coordinate extent of an image draw: (xmin, xmax, ymin, ymax)."
function image_extent(im::ImageRec)
    a, _, _, d, e, f = im.ctm
    return (min(e, e + a), max(e, e + a), min(f, f + d), max(f, f + d))
end

"Page coordinates of the pixel centres of an image draw: (x of each column, y of each row)."
function pixel_centres(im::ImageRec)
    a, b, c, d, e, f = im.ctm
    (abs(b) < 1e-9 && abs(c) < 1e-9) || error("rotated image placement not supported")
    xs = [e + a * (j - 0.5) / im.W for j in 1:im.W]
    ys = [f + d * (1 - (i - 0.5) / im.H) for i in 1:im.H]
    return xs, ys
end

# ---------------------------------------------------------------------------------------------
# Axes: frames, ticks, linear maps
# ---------------------------------------------------------------------------------------------

segs(p::PathRec) = [sp for sp in p.subpaths if length(sp) == 2]
isblack(c) = all(<(0.05), c)

"Black axis frame (the four 2-point spines) enclosing the point (x, y): (xmin, xmax, ymin, ymax)."
function frame_of(paths, x, y)
    vs = Tuple{Float64,Float64,Float64}[]
    for p in paths
        (p.kind == :stroke && isblack(p.color) && isempty(p.dash)) || continue
        for sp in segs(p)
            (x1, y1), (x2, y2) = sp
            if abs(x1 - x2) < 1e-3 && abs(y1 - y2) > 20     # vertical spine
                push!(vs, (x1, min(y1, y2), max(y1, y2)))
            end
        end
    end
    left = [v for v in vs if v[1] < x && v[2] < y < v[3]]
    right = [v for v in vs if v[1] > x && v[2] < y < v[3]]
    (isempty(left) || isempty(right)) && error("no frame around ($x, $y)")
    l = left[argmax(first.(left))]; r = right[argmin(first.(right))]
    return (l[1], r[1], max(l[2], r[2]), min(l[3], r[3]))
end

"short black vertical strokes hanging below y = ybot, inside xlims: tick x positions"
function ticks_below(paths, ybot, xlims; len = (0.3, 4.0))
    xs = Float64[]
    for p in paths, sp in segs(p)
        (p.kind == :stroke && isblack(p.color)) || continue
        (x1, y1), (x2, y2) = sp
        if abs(x1 - x2) < 1e-3 && len[1] < abs(y1 - y2) < len[2] && abs(max(y1, y2) - ybot) < 0.05 &&
           xlims[1] - 0.05 <= x1 <= xlims[2] + 0.05
            push!(xs, x1)
        end
    end
    return sort(unique(round.(xs; digits = 3)))
end

"short black horizontal strokes left of x = xleft, inside ylims: tick y positions"
function ticks_left(paths, xleft, ylims; len = (0.3, 4.0))
    ys = Float64[]
    for p in paths, sp in segs(p)
        (p.kind == :stroke && isblack(p.color)) || continue
        (x1, y1), (x2, y2) = sp
        if abs(y1 - y2) < 1e-3 && len[1] < abs(x1 - x2) < len[2] && abs(max(x1, x2) - xleft) < 0.05 &&
           ylims[1] - 0.05 <= y1 <= ylims[2] + 0.05
            push!(ys, y1)
        end
    end
    return sort(unique(round.(ys; digits = 3)))
end

"short black horizontal strokes starting at x ≈ xright (colourbar ticks), inside ylims"
function ticks_right(paths, xright, ylims; len = (0.3, 4.0), tol = 0.05)
    ys = Float64[]
    for p in paths, sp in segs(p)
        (p.kind == :stroke && isblack(p.color)) || continue
        (x1, y1), (x2, y2) = sp
        if abs(y1 - y2) < 1e-3 && len[1] < abs(x1 - x2) < len[2] && abs(min(x1, x2) - xright) < tol &&
           ylims[1] - 0.05 <= y1 <= ylims[2] + 0.05
            push!(ys, y1)
        end
    end
    return sort(unique(round.(ys; digits = 3)))
end

"""
    fit_axis(pos, vals) -> (f, maxres)

Least-squares linear map value = a + b·pos through the tick positions, and the worst tick residual
in value units (a measure of how exactly the ticks were placed).
"""
function fit_axis(pos::AbstractVector, vals::AbstractVector)
    length(pos) == length(vals) || error("$(length(pos)) tick positions for $(length(vals)) values")
    A = hcat(ones(length(pos)), pos)
    a, b = A \ vals
    res = maximum(abs.(a .+ b .* pos .- vals))
    return (x -> a + b * x), res
end

# ---------------------------------------------------------------------------------------------
# Colour-map inversion
# ---------------------------------------------------------------------------------------------

"""
    colorbar_lut(pix, im::ImageRec, fval) -> (lut, levels)

Colour → value table from a colourbar image (`pix`, H × W × 3) drawn by `im`, with `fval` the
fitted map page-y → value of its ticks. Each distinct colour is a contiguous band of rows; its value
is the mean of the band's row-centre values (the band centre). `levels` lists (colour, value) from
the bottom of the bar up.
"""
function colorbar_lut(pix, im::ImageRec, fval)
    _, ys = pixel_centres(im)
    col = size(pix, 2) ÷ 2
    bands = Dict{NTuple{3,UInt8},Vector{Float64}}()
    for i in axes(pix, 1)
        c = (pix[i, col, 1], pix[i, col, 2], pix[i, col, 3])
        push!(get!(bands, c, Float64[]), fval(ys[i]))
    end
    lut = Dict(c => sum(v) / length(v) for (c, v) in bands)
    levels = sort([(c, v) for (c, v) in lut]; by = last)
    return lut, levels
end

"""
    invert_image(pix, lut, levels) -> (vals, nmiss)

Value of every pixel through the colour table; colours missing from it (none expected) take the
nearest table colour in RGB, and are counted in `nmiss`.
"""
function invert_image(pix, lut, levels)
    H, W, _ = size(pix)
    vals = Matrix{Float64}(undef, H, W)
    nmiss = 0
    cols = [Float64.(collect(c)) for (c, _) in levels]
    for j in 1:W, i in 1:H
        c = (pix[i, j, 1], pix[i, j, 2], pix[i, j, 3])
        v = get(lut, c, nothing)
        if v === nothing
            nmiss += 1
            q = Float64[c...]
            k = argmin([sum(abs2, q .- cc) for cc in cols])
            v = levels[k][2]
        end
        vals[i, j] = v
    end
    return vals, nmiss
end

# ---------------------------------------------------------------------------------------------
# Lines
# ---------------------------------------------------------------------------------------------

"strokes of colour `rgb` (within tol) whose points all lie in the frame"
function strokes_in(paths, frame, rgb; tol = 0.01)
    xmin, xmax, ymin, ymax = frame
    out = PathRec[]
    for p in paths
        p.kind == :stroke || continue
        maximum(abs.(p.color .- rgb)) < tol || continue
        pts = reduce(vcat, p.subpaths; init = NTuple{2,Float64}[])
        isempty(pts) && continue
        all(xmin - 0.5 <= x <= xmax + 0.5 && ymin - 0.5 <= y <= ymax + 0.5 for (x, y) in pts) && push!(out, p)
    end
    return out
end

"all points of a stroke, in drawing order"
polyline(p::PathRec) = reduce(vcat, p.subpaths; init = NTuple{2,Float64}[])

end # module
