# Reader for the digitized SI maps of devices B and C (data_exp/digitized, refit_BC/digitize.jl).
#
# Returns the same `Data` as read_data (data_reader.jl), so prepare_data / extract_peak_currents /
# compute_normal_resistances / γAH run unchanged on it. Load after data_reader.jl.
#
# The one adaptation: the published maps are already R_s-corrected dV/dI, and their current axis is
# one image row per pixel (≈ 0.11 nA for device C's zoom, ≈ 0.72 nA for device B). Every smoothing
# and peak window of the pipeline is a number of *points*, tuned on the main device's run 275_16
# (step 0.2167 nA): Savitzky–Golay 50 points ≈ 10.8 nA in read_data, 101 points ≈ 21.9 nA in
# get_peak_indices_smoothed. The maps are therefore resampled (linearly) onto that same current step,
# so those windows keep their width in nA and no function of the pipeline needs a new parameter.

using DelimitedFiles

@with_kw struct Digitized_Data
    map::String                         # data_exp/digitized/devX_dVdI_*.dat
    line::String = ""                   # data_exp/digitized/devX_R0_line.dat (published R_0, checks only)
    δI::Float64 = 0.21666672            # nA, current step of run 250mK_DEV13_275_16 (main device, Fig. 3)
    Rs::Float64 = NaN                   # Ω, series resistance already subtracted in the published map
end

"digitized matrix file → (B (T), I (nA), M[iI, iB] dV/dI (Ω))"
function read_digitized_matrix(path)
    raw = readdlm(path, '\t', Float64; comments = true, comment_char = '#')
    return raw[1, 2:end], raw[2:end, 1], raw[2:end, 2:end]
end

"linear interpolation of y(x) (x sorted) at xq"
function lin_resample(x, y, xq)
    out = similar(xq, Float64)
    for (k, q) in enumerate(xq)
        j = clamp(searchsortedlast(x, q), 1, length(x) - 1)
        t = (q - x[j]) / (x[j+1] - x[j])
        out[k] = (1 - t) * y[j] + t * y[j+1]
    end
    return out
end

function read_data(d::Digitized_Data; exp_params = nothing, window_size = 50, order = 3)
    Bs, Iraw, M = read_digitized_matrix(d.map)
    p = sortperm(Iraw)
    Iraw, M = Iraw[p], M[p, :]
    # symmetric grid of step δI through I = 0, inside the digitized range
    n = floor(Int, min(-first(Iraw), last(Iraw)) / d.δI)
    Is = collect((-n:n) .* d.δI)
    dVdImat = Matrix{Float64}(undef, length(Bs), length(Is))
    for j in eachindex(Bs)
        dVdImat[j, :] .= lin_resample(Iraw, M[:, j], Is)
    end
    dVdImat = sg_filter_along_y(dVdImat; window_size, order)
    # Is in A, as read_data returns them (prepare_data rescales with Icoef = 1e9)
    return Data(; Ts = Bs, Is = Is .* 1e-9, dVdImat = dVdImat)
end

function prepare_data(d::Digitized_Data, Imax::Real = 120; Icoef = 1e9, window_size = 50, kw...)
    data = read_data(d; window_size)
    @unpack Ts, Is, dVdImat = data
    Is = Is .* Icoef
    idx = findall(Is .<= Imax .&& Is .>= -Imax)
    return Ts, Is[idx], dVdImat[:, idx], nothing, nothing
end

"published R_0(B) line of panel (c): (B, R0)"
function read_R0_line(d::Digitized_Data)
    raw = readdlm(d.line, '\t', Float64; comments = true, comment_char = '#')
    return raw[:, 1], raw[:, 2]
end
