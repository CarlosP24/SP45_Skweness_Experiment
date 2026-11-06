include("plotters/theo_exp_core.jl")

# Referee request: same junction as Fig. 3 (fitted_$(name)) but with a fully depleted
# core (depleted_$(name), only μ changed). CdGM analogs need occupied core subbands,
# so without them Ic just follows the shell's Little-Parks gap and loses its skewness.
function fig_depleted(name::String; R = 64, conv = 74.4, Blims = (0, 0.38))
    Φ0 = 2067.833 #T⋅nm^2
    Bns = [1, 2].* Φ0 / (π * R^2)
    Bhalfs = [0.5, 1.5].* Φ0 / (π * R^2)

    # B = 0 is dropped by calc_Josephson; the first point (≈ 2 mT) stands in for it
    Brng, Ic = ics("fitted_$(name)"; conv)
    Brng_d, Ic_d = ics("depleted_$(name)"; conv = 1e3 * conv)     # pA
    Ic0, Ic0_d = first(Ic), first(Ic_d)

    fig = Figure(size = (600, 330), fontsize = 20)

    ax = Axis(fig[1, 1]; xlabel = L"$B$ (T)", ylabel = L"$I_c / I_c(B = 0)$")
    lines!(ax, Brng, Ic ./ Ic0; color = :black, linewidth = 4)
    lines!(ax, Brng, Ic ./ Ic0; color = :orange, linewidth = 3,
        label = L"\text{Fitted to data}, \; I_c(0) = %$(round(Int, Ic0))\,\text{nA}")
    lines!(ax, Brng_d, Ic_d ./ Ic0_d; color = :black, linewidth = 4)
    lines!(ax, Brng_d, Ic_d ./ Ic0_d; color = :dodgerblue, linewidth = 3,
        label = L"\text{Depleted}, \; I_c(0) = %$(round(Int, Ic0_d))\,\text{pA}")
    xlims!(ax, Blims)
    ylims!(ax, 0, 1.15)
    ax.yticks = 0:0.5:1

    vlines!(ax, Bns; color = :black, linestyle = :dash)
    vlines!(ax, Bhalfs; color = :gray, linestyle = :dash)

    for n in 0:2
        text!(ax, (n + 0.25) * Φ0 / (π * R^2), 0.1; text = "$n"*"L", align = (:center, :center), color = :black)
    end

    axislegend(ax; position = :rt, labelsize = 16, rowgap = 0, padding = (6, 6, 2, 2))

    gax = Axis(fig[1, 1]; xaxisposition = :top, xlabel = L"\Phi / \Phi_0")
    xlims!(gax, Blims)
    gax.xticks = (vcat(Bns, Bhalfs), [L"1", L"2", L"0.5", L"1.5"])
    hidespines!(gax)
    hideydecorations!(gax)
    hidexdecorations!(gax, ticks = false, ticklabels = false, label = false)

    return fig
end
