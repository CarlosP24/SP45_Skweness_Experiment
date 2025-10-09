# helper: compute maximum while ignoring NaN values; return -Inf if all are NaN
safe_max(vals) = begin
    arr = collect(vals)
    arr = filter(x -> !isnan(x), arr)
    isempty(arr) ? -Inf : maximum(arr)
end

function fit_LDOS(name::String)
    wire_system = wires_fit[name]
    @unpack wire, params = wire_system
    @unpack µrng, τΓrng, wrng, grng, Δws, BsLDOS, ω, skew2, outdir = params 

    ΦtoB = get_B(wire)

    Brng1L = filter(B -> ΦtoB(0.501) <= B <= ΦtoB(1.0), BsLDOS)
    Brng1R = filter(B -> ΦtoB(1.0) < B <= ΦtoB(1.499), BsLDOS)

    Brng2L = filter(B -> ΦtoB(1.501) < B <= ΦtoB(2.0), BsLDOS)
    Brng2R = filter(B -> ΦtoB(2.0) < B <= ΦtoB(2.499), BsLDOS)

    if skew2 == "right"
        Brng2L, Brng2R = Brng2R, Brng2L
    end

    # Output path
    path = "$(outdir)/LDOS_fit/$(name).jld2"
    mkpath(dirname(path))


    # Build Hamiltonian
    hSM, hSC, wparams = build_cyl_mm(wire)

    # Build Green's function
    g_right, g = greens_dict["semi"](hSC, wparams)
    
    # Build LDOS
    ρ = ldos(g_right[cells = (-1,)])

    ρf(ω, B, µ, τΓ, w, ga) = try
        return sum(ρ(ω; B, ω, µ, τΓ, w = w, g = ga))
    catch
        return NaN
    end

    function ρparallel(µ, τΓ, w, Δw, ga)
        ρB = if Δw == 0
            B -> ρf(ω, B, µ, τΓ, w, ga)
        else
            B -> sum([ρf(ω, B, µ, τΓ, w + dw, ga) for dw in [0, Δw]])
        end
        
        # Criteria 1: skweness in first lobe (ignore NaNs)
        if safe_max(ρB.(Brng1L)) <= safe_max(ρB.(Brng1R))
            return false
        end

        # Criteria 2: skweness in second lobe (ignore NaNs)
        if safe_max(ρB.(Brng2L)) >= safe_max(ρB.(Brng2R))
            return false
        end

        return true
    end

    LDOSfits = pfunction(
        ρparallel,
        [µrng, τΓrng, wrng, Δws, grng]
    )

    pts = collect(Iterators.product(µrng, τΓrng, wrng, Δws, grng))
    true_indices = findall(LDOSfits)
    fitted_params = pts[true_indices]

    return Results(
        params = params,
        wire = wire,
        fitted_params = fitted_params,
        path = path
    )
end