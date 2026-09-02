# Ported from SP43 src/parallelizers/josephson.jl:31-48 (SHA dd7ad2f).
#
# NOTE: SP43's working tree had an uncommitted edit that emptied the `try` body
# of the *other* method in that file (`pjosephson(J, Brng, lg)`), leaving it
# returning `nothing`. This port is taken from the committed HEAD, not from that
# working tree, and SP43 was left untouched. See data/theory/PROVENANCE.md.

"""
    pcpr(f, Φrng, Zs, lg)

Evaluate `f(Φ, Z)` — which must return one `J_S(φ)` vector — over the product of
fluxes `Φrng` and conserved angular-momentum sectors `Zs`, returning
`Dict(Z => Matrix)` where each matrix is `length(Φrng) × length(φrng)`.

Kept un-summed on purpose: the sum over `Z` is the total supercurrent, and the
whole point of this calculation is the decomposition that summing destroys.

`lg` is the length of the phase range, needed so a failed point can return the
right-shaped `NaN` vector instead of poisoning the reshape.
"""
function pcpr(f, Φrng, Zs, lg::Int)
    pts = Iterators.product(Φrng, Zs)
    Jss = @showprogress pmap(pts) do pt
        Φ, Z = pt
        return try
            f(Φ, Z)
        catch e
            @warn "CPR failed at Φ=$Φ, Z=$Z. Output is NaN." exception = e
            fill(NaN, lg)
        end
    end
    arr = reshape(Jss, size(pts)...)
    # one (Φ × φ) matrix per sector
    return Dict([Z => mapreduce(permutedims, vcat, arr[:, i]) for (i, Z) in enumerate(Zs)])
end
