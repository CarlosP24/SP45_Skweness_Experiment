using NonlinearSolve
using LinearAlgebra
using Optimization
using OptimizationOptimJL
using NaNMath

ħ2ome = 76.1996
μBΦ0 = 119.6941183                  #Bohr magneton times magnetic flux quantum


function create_system(Φs, R, d, n, ms, g, max_mJ, min_mu, minΔ = 0.2)
    m = ms * ħ2ome  # effective mass
    RLP = R + d/2
    # Helper functions - make sure all types are compatible with AD
    ϕ(Rav, Φ) = n - Rav^2 / RLP^2 * Φ
    μf(mJ, μ, Rav, Φ) = μ - (4*mJ^2 + 1 + ϕ(Rav, Φ)^2) / (8 * m * Rav^2)  
    VZ(Rav, Φ) = 0.5 * ϕ(Rav, Φ) / (2 * m * Rav^2) + 0.5 * g * μBΦ0 * Φ / (π * RLP^2)
    CmJ(mJ, Rav) = -mJ / (2 * m * Rav^2)  
    AmJ(mJ, Rav, Φ) = -mJ * ϕ(Rav, Φ) / (2 * m * Rav^2)

    function f!(residual, u, p)
        # Use simple bounded transformations
        # u = [sigmoid(mJ/max_mJ), sigmoid(μ), logit((Rav - R/2)/(R/2)), log(Δ-minΔ)]
        
        mJ = max_mJ / (1.0 + exp(-u[1]))  # Ensures 0 < mJ < max_mJ
        μ = min_mu + exp(u[2])  # Ensures μ >= min_mu
        Rav = R/2 + (0.3*R) / (1.0 + exp(-u[3]))  # Ensures R/2 < Rav <= 0.8*R
        Δ = minΔ + exp(u[4])  # Ensures Δ >= minΔ
        
        sorted_Φs = sort(Φs, rev=true)
        for (i, Φ) in enumerate(sorted_Φs)
            ΔmJ = i - 1
            mJ_i = mJ - ΔmJ
            μ_term = μf(mJ_i, μ, Rav, Φ) - CmJ(mJ_i, Rav)
            A_term = AmJ(mJ_i, Rav, Φ) + VZ(Rav, Φ)
            residual[i] = μ_term^2 - A_term^2 + Δ^2
        end
        return residual
    end

    return f!
end

function round_mJ(mJ, n)
    """Round mJ to integer if n is odd, or to half-integer if n is even."""
    if isodd(n)
        return round(mJ)  # Round to nearest integer
    else
        return round(2*mJ) / 2  # Round to nearest half-integer
    end
end

function fit_tubular_parameters(Φs, R, d, n, m, g, minΔ = 0.2; u0=nothing, solver=nothing, verbose=false, use_optimization=true, max_mJ=20.0, min_mu=0.0)
    """
    Solve the nonlinear system to find tubular nanowire parameters.
    
    # Arguments
    - `Φs`: Array of flux values (must have at least 4 elements)
    - `R`: Inner radius parameter
    - `d`: Shell thickness
    - `n`: Winding number (integer or half-integer constraint on mJ)
    - `m`: Effective mass parameter
    - `g`: g-factor
    - `minΔ`: Minimum allowed gap parameter (default: 0.2)
    - `u0`: Initial guess in transformed space
    - `solver`: NonlinearSolve solver (default: tries multiple approaches)
    - `verbose`: Print solver diagnostics
    - `use_optimization`: Try optimization-based approach first
    - `max_mJ`: Maximum allowed value for angular momentum quantum number (default: 20.0)
    - `min_mu`: Minimum allowed value for chemical potential (default: 0.0)
    
    # Returns
    Named tuple with fields:
    - `mJ`: Angular momentum quantum number (rounded to integer/half-integer)
    - `μ`: Chemical potential
    - `Rav`: Average radius (positive, < R)
    - `Δ`: Gap parameter (>= 0.2)
    - `converged`: Whether the solver converged
    - `mJ_raw`: Unrounded mJ value from solver
    """
    
    # Smart initial guess based on realistic ranges:
    # 0 < mJ < max_mJ, μ >= min_mu, R/2 < Rav <= 0.8*R, Δ >= minΔ
    if isnothing(u0)
        # Start at middle of realistic ranges
        mJ_guess = min(7.0, max_mJ * 0.5)  # Half of max_mJ or 7, whichever is smaller
        μ_guess = max(20.0, min_mu + 20.0)  # min_mu + 20 or 20, whichever is larger
        Rav_guess = 0.65 * R  # 65% of R (between R/2 and 0.8*R)
        Δ_guess = max(1.0, minΔ + 0.8)  # minΔ + 0.8 or 1.0, whichever is larger
        
        # Transform to unconstrained space
        u0 = [
            log(mJ_guess / (max_mJ - mJ_guess)),  # logit-like transform for bounded mJ
            log(μ_guess - min_mu),  # log transform for μ >= min_mu
            0.0,  # logit-ish for Rav
            log(Δ_guess - minΔ)  # log transform for Δ >= minΔ
        ]
    end
    
    # Create the system
    f! = create_system(Φs, R, d, n, m, g, max_mJ, min_mu, minΔ)
    
    sol = nothing
    best_residual = Inf
    
    # Try optimization approach first (often more robust)
    if use_optimization
        try
            if verbose
                println("Trying least-squares optimization approach...")
            end
            
            # Objective: minimize sum of squared residuals
            function loss(u, p)
                resid = zeros(length(Φs))
                f!(resid, u, p)
                return sum(abs2, resid)
            end
            
            opt_func = OptimizationFunction(loss, Optimization.AutoForwardDiff())
            opt_prob = OptimizationProblem(opt_func, u0)
            opt_sol = solve(opt_prob, BFGS(); maxiters=10000)
            
            # Check the residual
            resid = zeros(length(Φs))
            f!(resid, opt_sol.u, nothing)
            resid_norm = norm(resid)
            
            if verbose
                println("  Optimization residual norm: ", resid_norm)
            end
            
            if resid_norm < best_residual
                best_residual = resid_norm
                sol = (u=opt_sol.u, resid=resid, retcode=(resid_norm < 1.0 ? :Success : :Default))
            end
        catch e
            if verbose
                println("  Optimization failed: ", e)
            end
        end
    end
    
    # Try NonlinearSolve solvers
    prob = NonlinearProblem(f!, isnothing(sol) ? u0 : sol.u)
    
    if isnothing(solver)
        solvers = [
            TrustRegion(; radius_update_scheme=RadiusUpdateSchemes.Bastin),
            LevenbergMarquardt(; damping_initial=100.0),
            LevenbergMarquardt(; damping_initial=1000.0),
            TrustRegion(; radius_update_scheme=RadiusUpdateSchemes.NLsolve),
            NewtonRaphson(),
        ]
        
        for s in solvers
            try
                if verbose
                    println("Trying solver: ", typeof(s))
                end
                temp_sol = solve(prob, s; abstol=1e-3, reltol=1e-3, maxiters=10000)
                resid_norm = norm(temp_sol.resid)
                
                if verbose
                    println("  Residual norm: ", resid_norm)
                    println("  Return code: ", temp_sol.retcode)
                end
                
                # Keep track of best solution
                if resid_norm < best_residual
                    sol = temp_sol
                    best_residual = resid_norm
                end
                
                if (temp_sol.retcode == :Success || temp_sol.retcode == :Default) && resid_norm < 1.0
                    if verbose
                        println("✓ Found approximation with ", typeof(s))
                    end
                    break
                end
            catch e
                if verbose
                    println("✗ Failed with ", typeof(s), ": ", e)
                end
                continue
            end
        end
    else
        sol = solve(prob, solver; abstol=1e-6, reltol=1e-6, maxiters=5000)
    end
    
    if isnothing(sol)
        error("All solvers failed completely. The system may not have a solution with these constraints.")
    end
    
    if best_residual > 5.0
        @warn "Poor convergence: residual norm = $best_residual. Results may not be reliable."
    elseif best_residual > 1.0
        @warn "Moderate convergence (approximation): residual norm = $best_residual"
    elseif best_residual < 0.1
        println("✓ Good approximation found: residual norm = $best_residual")
    end
    
    # Transform back from unconstrained to constrained space
    mJ_raw = max_mJ / (1.0 + exp(-sol.u[1]))
    mJ_rounded = round_mJ(mJ_raw, n)
    μ = min_mu + exp(sol.u[2])
    Rav = R/2 + (0.3*R) / (1.0 + exp(-sol.u[3]))
    Δ = minΔ + exp(sol.u[4])
    
    if verbose
        println("\n" * "="^50)
        println("Final residual norm: ", norm(sol.resid))
        println("Residual vector: ", sol.resid)
        println("Transformed variables: ", sol.u)
    end
    
    return (
        mJ = mJ_rounded,
        μ = μ,
        Rav = Rav,
        Δ = Δ,
        converged = best_residual < 1.0,
        mJ_raw = mJ_raw,
        residual = sol.resid,
        residual_norm = norm(sol.resid)
    )
end

## Example usage
Φ0 = 2067.83 # T·nm²
Bs = [0.295, 0.345, 0.39]
R = 60
d = 8
n = 2
m = 0.023
g = 17
τ = 2.5

Φs = Bs .* π .* (R + d/2)^2 / Φ0
# Try with verbose output to diagnose convergence issues
result = fit_tubular_parameters(Φs, R, d, n, m, g,  τ * 0.2; verbose=true, max_mJ=20, min_mu = 0)

# Display results
println("\n" * "="^50)
println("Solution:")
println("  mJ (raw)     = ", result.mJ_raw)
println("  mJ (rounded) = ", result.mJ, " (", isodd(n) ? "integer" : "half-integer", ")")
println("  μ            = ", result.μ)
println("  Rav          = ", result.Rav)
println("  Δ            = ", result.Δ)
println("\nConverged: ", result.converged)
println("Residual: ", result.residual)

## Fitting each Φ individually
function Φzeros(RLP, n, Γ, μ, m)
    Delta(mJ, Rav) = 4 * sqrt(2* m) * Rav * NaNMath.sqrt((1 - 2 * mJ)^2 * μ - 2 * m * Rav^2 * Γ^2)
    Alpha(mJ, Rav) = 1 + 4 * mJ * (mJ - 1) + 8 * m * Rav^2 * μ
    f(mJ, Rav) = (RLP / Rav)^2 * (n - NaNMath.sqrt(Alpha(mJ, Rav) - Delta(mJ, Rav)))
    return f
end 

Bs1 = [0.14, 0.17]
Φs1 = Bs1 .* π .* (R + d/2)^2 / Φ0
mJ = 7
f = Φzeros(R + d/2, n, 2.5 * 0.2, 25, 0.023/ħ2ome)
for Φ in sort(Φs1, rev=true) 
    Rs = find_zeros(Rav -> f(mJ, Rav) - Φ, R/2, R)
    println("Rs for Φ=$Φ, mJ = $mJ : ", Rs)
    mJ -= 1
end

Bs2 = [0.295, 0.345, 0.39]
Φs2 = Bs2 .* π .* (R + d/2)^2 / Φ0
mJ = 13/2
f = Φzeros(R + d/2, n, 2.5 * 0.2, 25, 0.023/ħ2ome)
for Φ in sort(Φs2, rev=true) 
    Rs = find_zeros(Rav -> f(mJ, Rav) - Φ, R/2, R)
    println("Rs for Φ=$Φ, mJ = $mJ : ", Rs)
    mJ -= 1
end