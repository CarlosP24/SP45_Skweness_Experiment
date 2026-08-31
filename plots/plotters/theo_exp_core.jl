function ics(name::String; basepath = "data/Js", conv = 74.4)
    # Load results
    path = joinpath(basepath, "$(name).jld2")
    res = load(path)["res"]

    @unpack Js, params = res
    @unpack Brng, φrng = params

    Ic = getindex(findmax(Js; dims = 2), 1) |> vec

    return Brng, Ic * conv
end

function process_data(exp_data::Exp_Data, name::String; T::Real = 460, Imax = 120, colorrange = (2.9, 3.2), Blims = (0.08, 0.4), pref = 0.8)
    # Process experimental data
    data = read_data(exp_data)
    Ts, Is, dVdImat, _, _ = prepare_data(exp_data, Imax)

    Isw, _, _ = extract_peak_currents(dVdImat, Is, Ts)
    RNs = compute_normal_resistances(dVdImat, Is, Isw)

    ilim = findmin(abs.(Ts .- Blims[1]))[2]
    imax = findmin(abs.(Ts .- Blims[2]))[2]

    i0 = findmin(abs.(Is))[2]

    rnorm = dVdImat ./ RNs
    r0 = rnorm[:, i0] |> vec
    γs = γAH.(abs.(r0))
    IcAH = 0.0209 * T .* γs  # T in mK

    Ic = vcat(NaN .* ones(1:ilim), IcAH[ilim+1:imax])

    Trng = Ts[1:imax]

    Rs = r0[ilim:imax]
    Trng_Rs = Ts[ilim:imax]

    # Process theory results
    Brng, Ict = ics(name)
    Ict .*= pref
    jmin = findmin(abs.(Brng .- Blims[1]))[2]
    jmax = findmin(abs.(Brng .- Blims[2]))[2]

    R(γ) = besseli0(γ / 2)^(-2)
    Rs_theo = R.(Ict[1:jmax] ./ (0.0209 * T))

    Brng_Rs = Brng[1:jmax]

    return Trng, Ic, Is, Brng, Ict, Rs, Trng_Rs, Rs, Brng_Rs, Rs_theo, rnorm[ilim:imax, :], Ic[ilim:imax], dVdImat[1:imax, :]
end

function panel_resistance(ax, Trng_Rs, Rs, Brng_Rs, Rs_theo; kw...)
    scatter!(ax, Trng_Rs, Rs; color = :black, label = "Experiment", markersize = 8, strokecolor = :white, strokewidth = 0.5, kw...)
    lines!(ax, Brng_Rs, Rs_theo; color = :black, linewidth = 4, kw...)
    lines!(ax, Brng_Rs, Rs_theo; color = :orchid, label = "Theory", linewidth = 3, kw...)
    axislegend(ax; position = (1, 0.5))
end

function panel_dVdI(ax, dVdImat, Ts, Is, Ic, Brng, Ict; theo_color=:orange, exp_color= :limegreen, legend = true, kw...)
    hmap = heatmap!(ax, Ts, Is, dVdImat; rasterize = 5, colormap = :davos, colorrange = (10^2.9, 10^3.2), colorscale = Makie.pseudolog10, kw...)
    scatter!(ax, Ts, Ic; color = exp_color, markersize = 8, label = L"I_c^\text{AH}", strokewidth = 0.5, strokecolor = :black)
    scatter!(ax, Ts, -Ic; color = exp_color, markersize = 8, strokecolor = :black, strokewidth = 0.5)

    lines!(ax, Brng, Ict; color = :black, linewidth = 4)
    lines!(ax, Brng, -Ict; color = :black, linewidth = 4)
    lines!(ax, Brng, Ict; color = theo_color, linewidth = 3, label = L"I_c^\text{th}")
    lines!(ax, Brng, -Ict; color = theo_color, linewidth = 3)
    lims = (first(Ts), last(Ts))
    xlims!(ax, lims)
    if legend
        axislegend(ax; position = :rt, )
    end

    return lims, hmap
end

function panel_LDOS(ax, name::String, side::String; basepath = "data/LDOS", colorrange = (0, 1e-1), kw...)
    # Load results
    path = joinpath(basepath, "fitted_$(name)_$(side).jld2")
    res = load(path)["res"]

    @unpack LDOS, params = res
    @unpack Brng, ωrng = params
    ωrng = real.(ωrng)
    #ωrng = vcat(ωrng, -reverse(ωrng)[2:end])
    #LDOS = cat(LDOSR, reverse(LDOSL, dims = 2)[:, 2:end], dims = 2)
    if side == "left"
        reverse!(LDOS, dims = 2)
        reverse!(ωrng)
        ωrng .*= -1
    end
    heatmap!(ax, Brng, ωrng, LDOS; colormap = :thermal, colorrange, rasterize = 5, kw...)
    hlines!(ax, [0]; color = :white, linewidth = 2)
end
