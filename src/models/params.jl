@with_kw struct Calc_Params
    Brng = range(0, 0.6, 401)
    ωrng = range(-0.3, 0, 401) .+ 1e-4im
    φrng = range(0, 2π, 101)
    outdir = "data"
end

@with_kw struct Results
    params = nothing
    wire = nothing
    junction = nothing
    LDOS = nothing
    Js = nothing
    path = nothing
end