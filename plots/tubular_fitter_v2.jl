using Optimization
using OptimizationOptimJL
using LinearAlgebra

ħ2ome = 76.1996
μBΦ0 = 119.6941183
πoΦ0 = 1.519267e-3 
R = 60
d = 8
RLP = R + d/2
m = 0.023 / ħ2ome
g = 17
Δ = 0.2
Bs = [0.11,  0.17, 0.17, 0.195, 0.24]

function f(B, RLP, m, g, Rav, μ, Γ, mJ)
    Φ = B * πoΦ0 * RLP^2
    n = round(Int, Φ)
    ϕ = n - (Rav^2 / RLP^2) * Φ
    μf = μ - (4 * mJ^2 + 1 + ϕ^2) / (8 * m * Rav^2)
    VZ = 0.5 * (ϕ / (2 * m * Rav^2) + g * μBΦ0 * Φ / (π * RLP^2))
    AmJ = -mJ * ϕ / (2 * m * Rav^2)
    CmJ = -mJ / (2 * m * Rav^2)
    return (μf - CmJ)^2 - (AmJ + VZ)^2 + Γ^2
end

"""
    get_n(B, RLP)

Compute winding number n for a given B field.
"""
get_n(B, RLP) = round(Int, B * πoΦ0 * RLP^2)

"""
    round_mJ(mJ_raw, n)

Round mJ to integer if n is odd, half-integer if n is even.
"""
function round_mJ(mJ_raw, n)
    if isodd(n)
        return round(mJ_raw)
    else
        return round(2 * mJ_raw) / 2
    end
end

"""
    fit_zeros(Bs, RLP, m, g, Δ, R; 
              max_mJ=20.0, B_tol=0.01, min_mu=0.0, max_attempts=10, target_residual=0.1, verbose=false)

Fit the zeros of function f to find optimal Rav, μ, Γ, individual mJs, and adjusted Bs.

# Constraints:
- 0.1*R < Rav < 0.8*R
- Γ > Δ
- μ >= min_mu
- mJ is integer if n is odd, half-integer if n is even
- Largest B has largest mJ, mJs decrease as B decreases (not necessarily consecutive)
- B values can vary within ±B_tol of input values, but n must remain fixed

# Multi-start:
- Tries up to max_attempts with different initial conditions until residual_norm < target_residual

# Returns:
Named tuple with fitted parameters and diagnostics.
"""
function fit_zeros(Bs, RLP, m, g, Δ, R; 
                   max_mJ=20.0, B_tol=0.01, min_mu=0.0, max_attempts=10, target_residual=0.1, verbose=false)
    
    nB = length(Bs)
    sorted_Bs_input = sort(Bs, rev=true)
    
    # Fix n values from input Bs (these are immutable)
    ns = [get_n(B, RLP) for B in sorted_Bs_input]
    
    # Loss function: sum of squared f values
    # Parameters: u[1]=Rav, u[2]=μ, u[3]=Γ, u[4:4+nB-1]=mJs, u[4+nB:end]=B_offsets
    function loss(u, p)
        # Parameter transformations to enforce constraints
        Rav = 0.1 * R + 0.7 * R / (1.0 + exp(-u[1]))  # 0.1*R < Rav < 0.8*R
        μ = min_mu + exp(u[2])               # μ >= min_mu
        Γ = Δ + exp(u[3])                    # Γ > Δ
        
        # mJs: each is bounded (0, max_mJ)
        mJs_raw = [max_mJ / (1.0 + exp(-u[3 + i])) for i in 1:nB]
        
        # B offsets: each is bounded (-B_tol, +B_tol) via tanh
        B_offsets = [B_tol * tanh(u[3 + nB + i]) for i in 1:nB]
        Bs_adj = sorted_Bs_input .+ B_offsets
        
        # Check for invalid values (NaN or Inf) and add huge penalty
        penalty = 0.0
        if any(isnan.(Bs_adj)) || any(isinf.(Bs_adj)) || any(isnan.(mJs_raw)) || any(isinf.(mJs_raw)) || 
           isnan(Rav) || isinf(Rav) || isnan(μ) || isinf(μ) || isnan(Γ) || isinf(Γ)
            return 1e10  # Return huge loss for invalid values
        end
        
        # Penalty for mJs_raw not decreasing with sufficient gap
        # Need at least 1.0 separation before rounding to ensure strict inequality after rounding
        for i in 2:nB
            gap = mJs_raw[i-1] - mJs_raw[i]
            if gap < 1.0  # Need at least 1 unit separation
                penalty += 10000.0 * (1.0 - gap)^2
            end
        end
        
        # Round mJs to integer/half-integer based on FIXED n parity
        mJs = [round_mJ(mJ_raw, n) for (mJ_raw, n) in zip(mJs_raw, ns)]
        
        # Additional penalty if rounded mJs are not strictly decreasing (should not happen with gap >= 1)
        for i in 2:nB
            if mJs[i] >= mJs[i-1]  # Equal or larger is not allowed
                penalty += 100000.0  # Huge penalty as this violates hard constraint
            end
        end
        
        # Compute sum of squared residuals
        total = 0.0
        for (B, mJ) in zip(Bs_adj, mJs)
            fval = f(B, RLP, m, g, Rav, μ, Γ, mJ)
            if isnan(fval) || isinf(fval)
                return 1e10  # Return huge loss if function evaluation fails
            end
            total += fval^2
        end
        
        return total + penalty
    end
    
    # Try multiple initial conditions until we get acceptable residual
    best_result = nothing
    best_residual_norm = Inf
    
    for attempt in 1:max_attempts
        if verbose
            println("\n" * "="^50)
            println("Attempt ", attempt, "/", max_attempts)
        end
        
        # Generate initial guess with some randomness after first attempt
        if attempt == 1
            # First attempt: use priors
            Rav_guess = 0.45 * R  # Middle of 0.1R to 0.8R
            μ_guess = max(20.0, min_mu + 20.0)
            Γ_guess = 2.5 * Δ
            mJ_guesses = [max_mJ - 2 * i for i in 0:(nB-1)]
        else
            # Subsequent attempts: add randomness
            Rav_guess = (0.1 + 0.7 * rand()) * R  # 0.1R to 0.8R
            μ_guess = min_mu + 20.0 * (0.5 + rand())  # min_mu+10 to min_mu+40
            Γ_guess = Δ * (2.0 + 2.0 * rand())     # 2Δ to 4Δ
            # Ensure mJ_guesses are strictly decreasing and in valid range (0, max_mJ)
            mJ_guesses = Float64[]
            for i in 0:(nB-1)
                mJ_guess = max_mJ * (0.4 + 0.4 * rand()) * (1.0 - 0.15 * i)
                # Clamp to valid range with margin
                mJ_guess = clamp(mJ_guess, 0.1, max_mJ - 0.1)
                push!(mJ_guesses, mJ_guess)
            end
        end
        
        # Ensure all mJ_guesses are valid for logit transform
        mJ_guesses = [clamp(mJ, 0.01, max_mJ - 0.01) for mJ in mJ_guesses]
        
        u0 = vcat(
            [log((Rav_guess - 0.1 * R) / (0.8 * R - Rav_guess))],  # logit for Rav in [0.1R, 0.8R]
            [log(μ_guess - min_mu)],                    # log for μ - min_mu
            [log(Γ_guess - Δ)],                         # log for Γ - Δ
            [log(mJ / (max_mJ - mJ)) for mJ in mJ_guesses],  # logit for each mJ
            randn(nB) * 0.1  # B offsets with small random perturbations
        )
        
        # Optimize with error handling
        opt_func = OptimizationFunction(loss, Optimization.AutoForwardDiff())
        opt_prob = OptimizationProblem(opt_func, u0)
        
        local opt_sol
        try
            opt_sol = solve(opt_prob, BFGS(); maxiters=10000)
        catch e
            if verbose
                println("Optimization failed: ", typeof(e), " - ", e.msg)
                println("Skipping this attempt...")
            end
            continue  # Skip to next attempt
        end
        
        # Transform back to physical parameters
        u = opt_sol.u
        Rav = 0.1 * R + 0.7 * R / (1.0 + exp(-u[1]))
        μ = min_mu + exp(u[2])
        Γ = Δ + exp(u[3])
        
        mJs_raw = [max_mJ / (1.0 + exp(-u[3 + i])) for i in 1:nB]
        B_offsets = [B_tol * tanh(u[3 + nB + i]) for i in 1:nB]
        Bs_fitted = sorted_Bs_input .+ B_offsets
        
        # Check for invalid values before proceeding
        if any(isnan.(Bs_fitted)) || any(isinf.(Bs_fitted)) || any(isnan.(mJs_raw)) || any(isinf.(mJs_raw)) ||
           isnan(Rav) || isinf(Rav) || isnan(μ) || isinf(μ) || isnan(Γ) || isinf(Γ)
            if verbose
                println("Warning: Invalid values encountered in iteration, skipping...")
            end
            continue  # Skip this attempt
        end
        
        # Use fixed n values (verify they haven't changed)
        ns_check = [get_n(B, RLP) for B in Bs_fitted]
        mJs = [round_mJ(mJ_raw, n) for (mJ_raw, n) in zip(mJs_raw, ns)]
        
        # Compute individual residuals
        residuals = [f(B, RLP, m, g, Rav, μ, Γ, mJ) for (B, mJ) in zip(Bs_fitted, mJs)]
        residual_norm = norm(residuals)
        
        if verbose
            println("Residual norm: ", residual_norm)
            if any(ns_check .!= ns)
                @warn "Warning: Some n values changed!" ns ns_check
            end
        end
        
        # Keep track of best result
        if residual_norm < best_residual_norm
            best_residual_norm = residual_norm
            best_result = (
                Rav = Rav,
                μ = μ,
                Γ = Γ,
                Bs_input = sorted_Bs_input,
                Bs_fitted = Bs_fitted,
                mJs = mJs,
                ns = ns,
                residuals = residuals,
                residual_norm = residual_norm,
                converged = residual_norm < target_residual,
                attempt = attempt
            )
        end
        
        # Check if we achieved target
        if residual_norm < target_residual
            if verbose
                println("✓ Target residual achieved!")
            end
            break
        end
    end
    
    # Print final results
    if verbose
        println("\n" * "="^50)
        println("Optimization completed after ", best_result.attempt, " attempt(s)")
        println("Best residual norm: ", best_result.residual_norm)
        println("\nFitted parameters:")
        println("  Rav = ", best_result.Rav)
        println("  μ   = ", best_result.μ)
        println("  Γ   = ", best_result.Γ)
        println("\nB and mJ values (B descending):")
        for (B_in, B_fit, mJ, n) in zip(best_result.Bs_input, best_result.Bs_fitted, best_result.mJs, best_result.ns)
            ΔB = B_fit - B_in
            println("  B = ", round(B_fit, digits=4), " (input: ", B_in, ", Δ = ", round(ΔB, digits=4), ") → n = ", n, ", mJ = ", mJ)
        end
        println("\nIndividual residuals:")
        for (B, r) in zip(best_result.Bs_fitted, best_result.residuals)
            println("  B = ", round(B, digits=4), " → f = ", r)
        end
        
        if !best_result.converged
            @warn "Did not achieve target residual after $(max_attempts) attempts. Best residual: $(best_result.residual_norm)"
        end
    end
    
    return best_result
end

# Example usage
result = fit_zeros(Bs, RLP, m, g, Δ, R; max_mJ=5.0, verbose=true, B_tol=0.02, max_attempts=1000, min_mu=30.0,)

println("\n" * "="^50)
println("Summary:")
println("  Rav = ", result.Rav)
println("  μ   = ", result.μ)
println("  Γ   = ", result.Γ)
println("  mJs = ", result.mJs)
println("  Converged: ", result.converged)