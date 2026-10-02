@with_kw struct Exp_Params
    Rs::Float64 = 5490
    β::Float64 = 0.95
    dI::Float64 = 5e-9 * β
    G_0::Float64 = 7.74809173e-5
end

@with_kw struct Exp_Data
    keysight = "data_exp/250mK_DEV13_275_16_dmm_keysight_signal.dat"
    lockin = "data_exp/250mK_DEV13_275_16_mf1_X.dat"
end

@with_kw struct Data
    Ts = Float64[]
    Is = Float64[]
    Vmat = Array{Float64,2}(undef, 0, 0)
    dVmat = Array{Float64,2}(undef, 0, 0)
    dIdVmat = Array{Float64,2}(undef, 0, 0)
    dVdImat = Array{Float64,2}(undef, 0, 0)
end

function sg_filter_along_y(mat::AbstractMatrix{<:Real}; window_size::Int = 9, order::Int = 3)
    out = copy(mat)
    nrows, ncols = size(out)

    # Ensure an odd window that fits the row length and polynomial order.
    w = min(window_size, ncols)
    w = iseven(w) ? w - 1 : w

    if w <= order || w < 3
        return out
    end

    for r in 1:nrows
        row = collect(Float64, @view out[r, :])
        if all(isfinite, row)
            out[r, :] .= SavitzkyGolay.savitzky_golay(row, w, order).y
        end
    end

    return out
end

function sg_filter_along_x(mat::AbstractMatrix{<:Real}; window_size::Int = 9, order::Int = 3)
    out = copy(mat)
    nrows, ncols = size(out)

    # Ensure an odd window that fits the column length and polynomial order.
    w = min(window_size, nrows)
    w = iseven(w) ? w - 1 : w

    if w <= order || w < 3
        return out
    end

    for c in 1:ncols
        col = collect(Float64, @view out[:, c])
        if all(isfinite, col)
            out[:, c] .= SavitzkyGolay.savitzky_golay(col, w, order).y
        end
    end

    return out
end

function read_data(exp_data::Exp_Data; exp_params::Exp_Params = Exp_Params(; ), window_size = 50, order = 3)
    @unpack keysight, lockin = exp_data
    @unpack Rs, dI, G_0 = exp_params

    # load files (DelimitedFiles.readdlm returns an Array)

    if keysight === nothing
        return read_data_nokeysight(exp_data; exp_params)
    end

    keysight = DelimitedFiles.readdlm(keysight)
    lockin = DelimitedFiles.readdlm(lockin)

    # unique sorted indices for rows/columns
    Ts = sort(unique(keysight[:, 1]))
    Is = sort(unique(keysight[:, 2]))

    # mapping from value -> matrix index
    t_idx = Dict(val => i for (i, val) in enumerate(Ts))
    i_idx = Dict(val => j for (j, val) in enumerate(Is))

    # create matrix with missing values for pairs not present
    Vmat = fill(NaN, length(Ts), length(Is))
    dVmat = fill(NaN, length(Ts), length(Is))
    dIdVmat = fill(NaN, length(Ts), length(Is))
    dVdImat = fill(NaN, length(Ts), length(Is))
    # populate matrix: column 1 -> row index, column 2 -> col index, column 4 -> value
    for r in eachindex(keysight[:, 1])
        t = keysight[r, 1]
        i = keysight[r, 2]
        v = keysight[r, 4]
        dv = abs.(lockin[r, 4])
        # Real voltage and lock-in correction
        Vmat[t_idx[t], i_idx[i]] = v - Is[i_idx[i]] * Rs
        dVmat[t_idx[t], i_idx[i]] = dv - dI * Rs
        dIdVmat[t_idx[t], i_idx[i]] = dI / dVmat[t_idx[t], i_idx[i]] / G_0
        dVdImat[t_idx[t], i_idx[i]] = dVmat[t_idx[t], i_idx[i]] / dI
    end

    dVdImat = sg_filter_along_y(dVdImat; window_size, order)

    return Data(; Ts = Ts, Is = Is, Vmat = Vmat, dVmat = dVmat, dIdVmat = dIdVmat, dVdImat = dVdImat)
end

function read_data_nokeysight(exp_data::Exp_Data; exp_params::Exp_Params = Exp_Params(;), window_size = 50, order = 3)
    @unpack lockin = exp_data
    @unpack Rs, dI, G_0 = exp_params

    lockin = DelimitedFiles.readdlm(lockin)

     # unique sorted indices for rows/columns
    Ts = sort(unique(lockin[:, 1]))
    Is = sort(unique(lockin[:, 2]))

    # mapping from value -> matrix index
    t_idx = Dict(val => i for (i, val) in enumerate(Ts))
    i_idx = Dict(val => j for (j, val) in enumerate(Is))

    # create matrix with missing values for pairs not present
    dVmat = fill(NaN, length(Ts), length(Is))
    dIdVmat = fill(NaN, length(Ts), length(Is))
    dVdImat = fill(NaN, length(Ts), length(Is))

    # populate matrix: column 1 -> row index, column 2 -> col index, column 4 -> value
    for r in eachindex(lockin[:, 1])
        t = lockin[r, 1]
        i = lockin[r, 2]
        dv = lockin[r, 4]

        # Real voltage and lock-in correction
        dVmat[t_idx[t], i_idx[i]] = dv - dI * Rs
        dIdVmat[t_idx[t], i_idx[i]] = dI / dVmat[t_idx[t], i_idx[i]] / G_0
        dVdImat[t_idx[t], i_idx[i]] = dVmat[t_idx[t], i_idx[i]] / dI
    end

    dVdImat = sg_filter_along_y(dVdImat; window_size, order)

    return Data(; Ts = Ts, Is = Is, dVmat = dVmat, dIdVmat = dIdVmat, dVdImat = dVdImat)
end


function prepare_data(exp_data::Exp_Data, Imax::Real=120; Icoef = 1e9, exp_params::Exp_Params = Exp_Params(; ), window_size = 50)
    """Load and preprocess experimental data."""
    data = read_data(exp_data; exp_params, window_size)
    @unpack Ts, Is, dVdImat, dVmat, Vmat = data 
    
    Is .*= Icoef
    idx = findall(Is .<= Imax .&& Is .>= -Imax) 
    Is = Is[idx]
    dVdImat = dVdImat[:, idx]
    dVmat = dVmat[:, idx]

    try
        Vmat = Vmat[:, idx]
    catch e 
    end
    
    return Ts, Is, dVdImat, dVmat, Vmat
end

