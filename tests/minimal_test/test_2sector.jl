# ═══════════════════════════════════════════════════════════════════════════════
# 2-sector diagnostic — understand the mobile labor model economics
# ═══════════════════════════════════════════════════════════════════════════════

using LinearAlgebra
using Printf

# ── 2-sector setup ──
N = 2
Ω = [0.5 0.5; 0.5 0.5]  # symmetric intermediates
factor_share = [0.5, 0.5]  # labor share
consumption_share = [0.5, 0.5]
supply_shock = [1.0, 1.0]
demand_shock = [1.5, 1.0]  # 50% demand increase in sector 1
labor_bar = 1.0  # total baseline labor

θ = 0.5  # elasticity between goods
ϵ = 0.5  # elasticity between labor and intermediates
σ = 0.9  # consumption elasticity

println("="^70)
println("  2-Sector Mobile Labor Model")
println("="^70)
println("  Ω = $Ω")
println("  factor_share = $factor_share")
println("  consumption_share = $consumption_share")
println("  demand_shock = $demand_shock")
println("  θ=$θ, ϵ=$ϵ, σ=$σ")

# ── The problem function ──
function problem_2sec(X, params)
    (; Ω, consumption_share, factor_share, supply_shock, demand_shock, labor_bar, θ, ϵ, σ, η) = params
    N = 2

    p = max.(X[1:N], 1e-10)
    y = max.(X[N+1:2N], 1e-10)
    w = max(X[2N+1], 1e-10)

    out = zeros(eltype(X), 2N + 1)

    # Intermediate price
    IP = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

    # CPI
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Labor demand
    L_i = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ

    # Final demand
    total_income = w * sum(L_i)
    final_demand = (total_income * p .^ (-σ) .* demand_shock .* consumption_share) ./ cpi .^ (-σ)

    # Intermediary demand
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* IP .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # Cost function
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* IP .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    # Eq 1: Price equation for sector 2 only (sector 1 = numeraire via CPI)
    out[1] = p[2] - cost[2]

    # Eq 2: Market clearing
    out[2:3] = y .- intermediary_demand .- final_demand

    # Eq 3: Labor market clearing
    out[4] = sum(L_i) - labor_bar * w^η

    # Eq 4: Numeraire (CPI = 1)
    out[5] = cpi - 1.0

    return out
end

# ── Newton solver (same as before) ──
function numerical_jacobian(f, x, params; eps=1e-6)
    n = length(x)
    f0 = f(x, params)
    m = length(f0)
    J = zeros(m, n)
    for j in 1:n
        x_plus = copy(x)
        x_plus[j] += eps
        f_plus = f(x_plus, params)
        J[:, j] = (f_plus - f0) / eps
    end
    return J
end

function solve_newton(f, x0, params; max_iter=500, tol=1e-10, verbose=false)
    x = copy(x0)
    n = length(x)

    for iter in 1:max_iter
        F = f(x, params)
        norm_F = norm(F)

        if verbose && (iter % 20 == 0 || iter == 1)
            @printf("    iter %d: ||F|| = %.2e\n", iter, norm_F)
        end

        if norm_F < tol
            return x, true, iter
        end

        J = numerical_jacobian(f, x, params)

        try
            dx = -(J \ F)
            alpha = min(1.0, 1.0 / max(1.0, norm(dx) / norm(x)))
            x_new = x + alpha * dx
            F_new = f(x_new, params)
            if norm(F_new) < norm_F
                x = x_new
            else
                x = x + 0.01 * alpha * dx
            end
        catch e
            x .+= 1e-4 * randn(n)
        end
    end

    return x, false, max_iter
end

# ── Run sweep ──
println("\n  η        p1       p2       y1       y2       w        L_total  GDP      Conv")
println("  " * "-"^85)

base_params = (
    Ω = Ω, consumption_share = consumption_share, factor_share = factor_share,
    supply_shock = supply_shock, demand_shock = demand_shock, labor_bar = labor_bar,
    θ = θ, ϵ = ϵ, σ = σ,
)

η_values = [0.0, 0.01, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
x0 = [1.0, 1.0, 0.5, 0.5, 1.0]  # p1, p2, y1, y2, w

for η in η_values
    params = merge(base_params, (η = η,))
    x_sol, converged, iters = solve_newton(problem_2sec, x0, params, verbose=false)

    p1, p2, y1, y2, w = x_sol

    # Compute derived quantities
    L_i = (max.(x_sol[1:2], 1e-10) .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (max.(x_sol[3:4], 1e-10) .^ (1 / ϵ)) ./ w) .^ ϵ
    L_total = sum(L_i)
    cpi = sum(consumption_share .* max.(x_sol[1:2], 1e-10) .^ (1 - σ))^(1 / (1 - σ))
    total_income = w * L_total
    cons = total_income .* consumption_share .* demand_shock .* (max.(x_sol[1:2], 1e-10) / cpi) .^ (-σ)
    gdp = sum(cons) / sum(consumption_share)

    @printf("  %-8.2f %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f  %s\n",
        η, p1, p2, y1, y2, w, L_total, gdp, converged ? "✅" : "⚠️")

    global x0 = copy(x_sol)
end

println("\n" * "="^70)
println("  Key insight: does w respond to η?")
println("="^70)