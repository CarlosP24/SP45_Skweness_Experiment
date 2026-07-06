function calc_LDOS(name::String)
    wire_system = wire_systems[name]
    @unpack wire, params, Δw = wire_system
    @unpack Brng, ωrng, outdir = params

    # Output path
    path = "$(outdir)/LDOS/$(name).jld2"
    mkpath(dirname(path))

    # Build Hamiltonian
    hSM, hSC, wparams = build_cyl_mm(wire)

    # Build Green's function
    g_right, g = greens_dict["semi"](hSC, wparams)

    ws = [wire.w]
    if Δw != 0
        push!(ws, wire.w + Δw)
    end

    # Calculate LDOS
    ρ = ldos(g_right[cells = (-1,)])
    LDOS = pfunction(
        (B, ω, w) -> try
            sum(ρ(ω; B, ω, w))
        catch e
            @warn "Error computing LDOS at (B=$B, ω=$ω): $e"
            NaN
        end,
        [Brng, ωrng, ws]
    )

    LDOS_sum = dropdims(sum(LDOS, dims = 3), dims = 3)

    return Results(
        params = params,
        wire = wire,
        LDOS = LDOS_sum,
        path = path
    )
end