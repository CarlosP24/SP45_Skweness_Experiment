# Check of the single-sector route used for devices B and C, on the main device (device A):
#  (a) find_τ_sm (Σ_Z single-sector normal conductance) against the multimode find_τ of fitted_C3;
#  (b) Σ_Z I_c at the published Zeeman convention (g_1D = 17π) against data/Js/fitted_C3.jld2;
#  (c) the same fluxes at the physical g = 17, and timing per (B, Z) point.
# Run from the SP45 root:  SP45_WORKERS=8 julia --project=. refit_BC/validate.jl
# Output: refit_BC/out/validate.txt
include(joinpath(pwd(), "refit_BC", "load.jl"))
using Printf

hero = junctions["fitted_C3"]
τ_mm = find_τ(hero)
const Zs = -9:9
pars = (μ = hero.wireL.μ, w = hero.wireL.w, τΓ = hero.wireL.τΓ)
wL17π, wR17π = device_wires("A"; pars..., g = 17π)
wL17, wR17 = device_wires("A"; pars..., g = 17.0)
τ_sm = find_τ_sm(wL17, wR17, hero.TN; Zs)

Bs = [0.04, 0.12, 0.16, 0.2, 0.28, 0.32]
chunks = [[b] for b in Bs]
t17π = @elapsed r17π = pmap(c -> js_reduced(wL17π, wR17π, τ_mm, c; T = hero.T, Zs), chunks)
t17 = @elapsed r17 = pmap(c -> js_reduced(wL17, wR17, τ_mm, c; T = hero.T, Zs), chunks)
Ic17π = [maximum(r.Jtot[1, :]) for r in r17π]
Ic17 = [maximum(r.Jtot[1, :]) for r in r17]

h = load("data/Js/fitted_C3.jld2")["res"]
Bh = collect(h.params.Brng)
Ich = vec(getindex(findmax(h.Js; dims = 2), 1))
ok = findall(isfinite, Ich)
Ic_hero = linear_interpolation(Bh[ok], Ich[ok]; extrapolation_bc = Line()).(Bs)

mkpath(joinpath("refit_BC", "out"))
open(joinpath("refit_BC", "out", "validate.txt"), "w") do io
    @printf(io, "τ: multimode find_τ %.5f, single-sector Σ_Z find_τ_sm %.5f (T_N = %.2f)\n", τ_mm, τ_sm, hero.TN)
    @printf(io, "timing: %.2f s per (B, Z) point per worker at g = 17π, %.2f s at g = 17 (%d workers)\n",
        t17π * nworkers() / (length(Bs) * length(Zs)), t17 * nworkers() / (length(Bs) * length(Zs)), nworkers())
    println(io, "B (T)   Ic_hero(mm)  Ic_sm(17π)  dev/peak   Ic_sm(17 phys)  ratio")
    for i in eachindex(Bs)
        @printf(io, "%-7.3f %-12.5f %-11.5f %+-10.2e %-15.5f %.3f\n", Bs[i], Ic_hero[i], Ic17π[i],
            (Ic17π[i] - Ic_hero[i]) / maximum(Ich[ok]), Ic17[i], Ic17[i] / Ic17π[i])
    end
end
print(read(joinpath("refit_BC", "out", "validate.txt"), String))
