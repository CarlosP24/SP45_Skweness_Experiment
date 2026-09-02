# m_J-resolved current-phase relation.
#
# Structurally this is calculations/Josephson.jl with three changes:
#
#   * wires come from `build_cyl` (m_J is the parameter `Z`) instead of
#     `build_cyl_mm` (m_J is a lattice coordinate);
#   * the sweep is over (Φ, Z) rather than (B, w), and the flux replaces the
#     field as the model's magnetic parameter;
#   * the result is kept as a Dict over Z instead of being traced away.
#
# The parameter conversions live in models/cpr.jl. Read that header before
# touching anything here.

"""
    calc_CPR(name::String)

Compute `J_S(φ)` resolved by conserved angular momentum `Z`, for the CPR system
`name`, and return the un-summed `Dict(Z => Φ × φ matrix)`.
"""
function calc_CPR(name::String)
    system = cpr_systems[name]
    @unpack wireL, wireR, Zs, Φrng, φrng, j_params, TN, kB, T, source, outdir = system
    @unpack imshift, atol, maxevals, order = j_params

    # Delete pathological points, exactly as the hero calculation does
    Φrng = collect(Φrng)
    deleteat!(Φrng, findall(x -> isapprox(x, 0), Φrng))
    φrng = filter_πs(φrng)

    path = "$(outdir)/$(name).jld2"
    mkpath(dirname(path))

    # Build the two single-mode wires and the junction Green's function
    hSM_left,  hSC_left,  params_left  = build_cyl(wireL)
    hSM_right, hSC_right, params_right = build_cyl(wireR)
    g_right, g_left, g = greens_dict["semi"](hSC_left, hSC_right, params_left, params_right;)

    # Find τ from TN. Recalibrated in this model rather than carried over: TN is
    # the fitted *normal transmission*, and τ is only the knob that reaches it,
    # so the transferable quantity is TN.
    τrng = subdiv(0, 1, 100)
    Gτs = get_TN(hSM_left, hSM_right, params_left, params_right, "semi", τrng; Φ = 0, Z = 0)
    Gτs = Gτs ./ maximum(Gτs)
    Tτ = linear_interpolation(τrng, Gτs)
    τ = find_zeros(τ -> Tτ(τ) - TN, 0, 1) |> first
    println("Found τ = $τ for TN = $TN")

    # Integration path — a function of Φ here, where the hero's is a function of B
    bw = maximum([bandwidth(params_left), bandwidth(params_right)])
    itipL = get_itip(params_left)
    itipR = get_itip(params_right)
    itip(x) = minimum([itipL(x), itipR(x)])
    ΩL = get_Ω(params_left)
    ΩR = get_Ω(params_right)
    Ω(x) = maximum([ΩL(x), ΩR(x)])
    Ωmin(x) = minimum([ΩL(x), ΩR(x)])
    ipath = Paths.polygon((mu, kBT; Φ = 0, _...) -> (-bw, -Ω(Φ), -Ω(Φ)/2 + itip(Φ)*1im, 0) .+ imshift*1im)

    J = josephson(g[attach_link["semi"]], ipath; omegamap = ω -> (; ω),
                  phases = φrng, atol, maxevals, order)

    # The gap closes at the lobe edges; the integrator has nothing to integrate
    # there. Same guard as the hero, in flux rather than field.
    function Jsector(Φ, Z)
        Ωmin(Φ) < 1e-3 && return fill(NaN, length(φrng))
        return J(kB * T; Φ, Z, τ)
    end

    Js = pcpr(Jsector, Φrng, Zs, length(φrng))

    write_cpr_sidecar(path, system, τ)

    return Results(
        params = system,
        junction = junctions[source],
        Js = Js,
        path = path
    )
end

# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

jsonesc(s) = replace(string(s), '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n")
jsonval(x::Bool) = string(x)
jsonval(x::Real) = isfinite(x) ? string(x) : "null"
jsonval(x::Complex) = jsonval(real(x))          # every Δ0 here is real
jsonval(x::Union{AbstractVector,AbstractRange}) = "[" * join(jsonval.(collect(x)), ", ") * "]"
jsonval(x) = "\"$(jsonesc(x))\""                # strings, symbols, everything else

"Read-only peek at SP43, for the record. Never writes to it."
function sp43_state()
    root = get(ENV, "SP43_PATH", joinpath(dirname(pwd()), "SP43_SNS_Paper"))
    isdir(joinpath(root, ".git")) || return (path = root, sha = "not found", dirty = "unknown")
    sha = try
        readchomp(`git -C $root rev-parse HEAD`)
    catch
        "unavailable"
    end
    dirty = try
        isempty(readchomp(`git -C $root status --porcelain`)) ? "clean" : "DIRTY"
    catch
        "unavailable"
    end
    return (path = root, sha = sha, dirty = dirty)
end

"""
    write_cpr_sidecar(path, system, τ)

Write the JSON sidecar §2 of the build doc asks for, next to the `.jld2`: the
full parameter set actually used, the SP43 commit the machinery came from and
whether that tree was dirty, the calibrated τ, and the date.

Hand-rolled rather than pulling in a JSON dependency for one flat object.
"""
function write_cpr_sidecar(path::AbstractString, system, τ)
    sp43 = sp43_state()
    wires = Dict("wireL" => system.wireL, "wireR" => system.wireR)
    fields = (:R, :w, :d, :a0, :t, :μ, :τΓ, :g, :Δ0, :ξd, :θ, :α,
              :Vmax, :Vmin, :ishollow, :shell, :L)

    io = IOBuffer()
    println(io, "{")
    println(io, "  \"generated\": ", jsonval(string(Dates_now())), ",")
    println(io, "  \"source_junction\": ", jsonval(system.source), ",")
    println(io, "  \"julia\": ", jsonval(string(VERSION)), ",")
    println(io, "  \"sp43\": {")
    println(io, "    \"path\": ", jsonval(sp43.path), ",")
    println(io, "    \"sha\": ", jsonval(sp43.sha), ",")
    println(io, "    \"worktree\": ", jsonval(sp43.dirty))
    println(io, "  },")
    println(io, "  \"calibrated_tau\": ", jsonval(τ), ",")
    println(io, "  \"TN\": ", jsonval(system.TN), ",")
    println(io, "  \"T_mK\": ", jsonval(system.T), ",")
    println(io, "  \"kB_meV_per_mK\": ", jsonval(system.kB), ",")
    println(io, "  \"Zs\": ", jsonval(collect(system.Zs)), ",")
    println(io, "  \"n_phi\": ", jsonval(length(system.Φrng)), ",")
    println(io, "  \"phi_min\": ", jsonval(minimum(system.Φrng)), ",")
    println(io, "  \"phi_max\": ", jsonval(maximum(system.Φrng)), ",")
    println(io, "  \"zeeman_convention\": ", jsonval(
        "g stored here is g_1D = pi * g_multimode; see src/models/cpr.jl"), ",")
    for (name, w) in sort(collect(wires), by = first)
        println(io, "  \"$name\": {")
        entries = ["    \"$f\": $(jsonval(getfield(w, f)))" for f in fields]
        println(io, join(entries, ",\n"))
        println(io, "  },")
    end
    println(io, "  \"note\": ", jsonval(
        "J_S(phi) is stored un-summed as Dict(Z => Phi x phi matrix) in the .jld2"))
    println(io, "}")

    open(splitext(path)[1] * ".json", "w") do f
        write(f, String(take!(io)))
    end
    return nothing
end

# Stdlib Dates is not a declared dependency of this project; this keeps the
# sidecar timestamped without adding one.
Dates_now() = Libc.strftime("%Y-%m-%dT%H:%M:%S%z", time())

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

"""
    verify_against_hero(res; rtol = 0.05)

Check that summing the m_J-resolved critical current reproduces the hero's, which
is the claim the whole route rests on: that the single-mode model with the
converted parameters is the same physical system as the fitted multimode one.

Returns a NamedTuple with the two Ic curves on a common flux axis and their
relative deviation. A disagreement here means a convention was missed — look at
the Zeeman factor of π first — and the figure must not be built until it passes.

This runs at T = 0 while the hero's cache is at 460 mK, so a smooth suppression
of the single-mode Ic at the few-percent level is expected and not a failure; a
*shape* difference, or a mismatch that grows across the lobe, is.
"""
function verify_against_hero(res; hero_path = "data/Js/fitted_C3.jld2", rtol = 0.05)
    system = res.params
    Φrng = collect(system.Φrng)
    deleteat!(Φrng, findall(x -> isapprox(x, 0), Φrng))

    # Σ_Z Ic^(Z): sum the sectors, then maximise over φ. Summing first matters —
    # the sectors peak at different phases, so max-then-sum would overcount.
    Jtot = sum(values(res.Js))
    Ic_1D = vec(maximum(Jtot; dims = 2))

    hero = load(hero_path)["res"]
    Ic_mm_all = vec(getindex(findmax(hero.Js; dims = 2), 1))
    Φ_mm = flux_of_B.(Ref(res.junction.wireL), collect(hero.params.Brng))

    # Restrict the hero to this calculation's flux window and interpolate onto it
    inwindow = findall(Φ -> first(Φrng) <= Φ <= last(Φrng), Φ_mm)
    itp = linear_interpolation(Φ_mm[inwindow], Ic_mm_all[inwindow]; extrapolation_bc = Line())
    Ic_mm = itp.(Φrng)

    scale = maximum(skipmissing(filter(!isnan, Ic_mm)))
    dev = abs.(Ic_1D .- Ic_mm) ./ scale
    worst = maximum(filter(!isnan, dev))

    ok = worst <= rtol
    ok ? @info("Σ_Z Ic matches the hero to $(round(100worst; digits = 2)) % of peak") :
         @warn("Σ_Z Ic deviates from the hero by up to $(round(100worst; digits = 2)) % of peak — " *
               "check the Zeeman π conversion in models/cpr.jl before building the figure")

    return (; Φ = Φrng, Ic_1D, Ic_mm, deviation = dev, worst, ok)
end
