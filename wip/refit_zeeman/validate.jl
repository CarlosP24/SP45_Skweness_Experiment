# Step 1 of the refit: does Σ_Z (single mode) reproduce the published multimode I_c?
# (a) single mode at g_1D = 17π (the published convention) vs data/Js/fitted_C3.jld2 at a few
#     fluxes in lobes 0-2 -> validates the sector route and times one point;
# (b) the same fluxes at the physical g = 17 -> first look at what the correction changes.
# Output: wip/refit_zeeman/validate.{jld2,txt}
include(joinpath(pwd(), "wip", "refit_zeeman", "load_sp45.jl"))

hero = junctions["fitted_C3"]
τ = find_τ(hero)                         # multimode, B = 0: unaffected by the Zeeman bug
@info "τ from the hero's multimode calibration" τ

Φs = [0.25, 0.75, 1.0, 1.25, 1.75, 2.0]
wL_pub, wR_pub = single_mode(hero.wireL), single_mode(hero.wireR)        # g_1D = 17π
wL_phy, wR_phy = sm_wire(hero.wireL; g = 17), sm_wire(hero.wireR; g = 17)  # g_1D = 17

t_pub = @elapsed res_pub = ic_single_mode(wL_pub, wR_pub, τ, Φs; T = hero.T)
t_phy = @elapsed res_phy = ic_single_mode(wL_phy, wR_phy, τ, Φs; T = hero.T)

# hero on the same fluxes
h = load("data/Js/fitted_C3.jld2")["res"]
Bh = collect(h.params.Brng)
Ich = vec(getindex(findmax(h.Js; dims = 2), 1))
Φh = flux_of_B.(Ref(hero.wireL), Bh)
ok = findall(isfinite, Ich)
itp = linear_interpolation(Φh[ok], Ich[ok]; extrapolation_bc = Line())
Ic_hero = itp.(Φs)

jldsave("wip/refit_zeeman/validate.jld2"; Φs, τ, Ic_hero, Ic_pub = res_pub.Ic, Ic_phy = res_phy.Ic,
    t_pub, t_phy, nw = nworkers())
open("wip/refit_zeeman/validate.txt", "w") do io
    println(io, "τ = $τ, workers = $(nworkers()), wall time: published g $(round(t_pub; digits=1)) s, physical g $(round(t_phy; digits=1)) s")
    println(io, "Φ/Φ0    Ic_hero(mm)   Ic_sm(g1D=17π)  rel.dev   Ic_sm(g=17 phys)  ratio phys/pub")
    for i in eachindex(Φs)
        println(io, rpad(Φs[i], 8), rpad(round(Ic_hero[i]; sigdigits=5), 14), rpad(round(res_pub.Ic[i]; sigdigits=5), 16),
            rpad(round((res_pub.Ic[i] - Ic_hero[i]) / maximum(Ic_hero); sigdigits=3), 10),
            rpad(round(res_phy.Ic[i]; sigdigits=5), 18), round(res_phy.Ic[i] / res_pub.Ic[i]; sigdigits=4))
    end
end
print(read("wip/refit_zeeman/validate.txt", String))
