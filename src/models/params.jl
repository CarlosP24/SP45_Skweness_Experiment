@with_kw struct Calc_Params
    Brng = range(0, 0.6, 401)
    ωrng = range(-0.3, 0, 401) .+ 1e-4im
    φrng = range(0, 2π, 101)
    outdir = "data"
end

@with_kw struct Fit_Params
    µrng = range(0, 50, step = 1)
    τΓrng = range(1, 5, step = 0.5)
    wrng = range(5, 30, step = 5)
    grng = range(0, 20, step = 5)
    φrng = range(0, 2π, 11)
    Δws = range(0, 30, step = 5)
    Bs01 = [0.01, 0.1,  0.2]       # B points to evaluate fitting conditions [0th lobe, 1st lobe]
    Bs1 = [0.1, 0.2]
    Bs2 = [0.28, 0.38]
    Idrop01 = 0.5            # Relative drop in I_c 0th to 1st lobe for fitting conditions
    Idrop01l = 0.05         # Relative drop in I_c 0th to 1st lobe, left side only
    Ibump1 = 2.0            # Relative drop in I_c skewed 1st lobe
    skew2 = "left"          # Second lobe skewness condition for fitting ("left", "right")
    outdir = "data"
    BsLDOS = range(0, 0.4, 101) 
    ω = 0.0 + 1e-3im
    prevfit = "data/LDOS_fit/C_left.jld2"  # Previous fit to start from
end

@with_kw struct Results
    params = nothing
    wire = nothing
    junction = nothing
    LDOS = nothing
    Js = nothing
    fitted_params = nothing
    path = nothing
end