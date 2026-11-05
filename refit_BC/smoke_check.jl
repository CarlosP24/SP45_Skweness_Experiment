# Sector-range check of the smoke runs: share of I_c carried by the sectors outside the range the full
# scan uses (B: |Z| ≤ 9, C: |Z| ≤ 12), at the largest μ and R_av of each grid.
#   julia --project=plots refit_BC/smoke_check.jl
using JLD2, Printf
for (name, zkeep) in (("smokeB", 9), ("smokeC", 12))
    f = joinpath("refit_BC", "out", "$(name)_es.jld2")
    isfile(f) || (println("$f missing"); continue)
    d = load(f)
    Zs, Bs = d["Zs"], d["Bs"]
    @printf("%s: candidate %s, τ = %.4f, wall %.2f h on %d workers\n", name, string(d["cands"][1]), d["τs"][1], d["wall_hours"], d["nworkers"])
    for (i, B) in enumerate(Bs)
        Ic = d["Ic"][1, i]
        JZ = d["JZc"][1, i, :]
        out = [abs(JZ[k]) for k in eachindex(Zs) if abs(Zs[k]) > zkeep]
        nz = [Zs[k] for k in eachindex(Zs) if abs(JZ[k]) > 1e-3 * abs(Ic)]
        @printf("  B = %.3f T: I_c = %.4f, outside |Z| ≤ %d: %.2e of I_c; sectors above 1e-3 I_c: %s\n",
            B, Ic, zkeep, sum(out) / abs(Ic), isempty(nz) ? "none" : "$(minimum(nz))…$(maximum(nz))")
    end
end
