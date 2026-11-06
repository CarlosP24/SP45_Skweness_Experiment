"""
    calc_LDOS(name)

LDOS at the end cell of the semi-infinite lead `name`, summed over orbitals and over the m_J sectors
`Zs` of the single-sector wire (`build_cyl`), on the grid Brng × ωrng (flux Φ = πoΦ0 B R_LP²).
"""
function calc_LDOS(name::String)
    wire_system = wire_systems[name]
    @unpack wire, params, Δw = wire_system
    @unpack Brng, ωrng, outdir = params

    # Output path
    path = "$(outdir)/LDOS/$(name).jld2"
    mkpath(dirname(path))

    Φ = get_Φ(wire)
    Zs = collect(wire.Zs)

    # The core width is a lattice property of build_cyl, so each w is its own build
    ws = [wire.w]
    if Δw != 0
        push!(ws, wire.w + Δw)
    end

    LDOS_sum = zeros(length(Brng), length(ωrng))
    for w in ws
        # Build Hamiltonian and Green's function
        hSM, hSC, wparams = build_cyl(Params(wire; w))
        g_right, g = greens_dict["semi"](hSC, wparams)
        ρ = ldos(g_right[cells = (-1,)])

        # One field per task: all energies and sectors of a field are cheap and done serially
        LDOS = pfunction(
            B -> [sum(Z -> try
                    sum(ρ(ω; Φ = Φ(B), ω, Z))
                catch e
                    @warn "Error computing LDOS at (B=$B, ω=$ω, Z=$Z): $e"
                    NaN
                end, Zs) for ω in ωrng],
            [Brng]
        )
        LDOS_sum .+= mapreduce(permutedims, vcat, LDOS)
    end

    return Results(
        params = params,
        wire = wire,
        LDOS = LDOS_sum,
        path = path
    )
end
