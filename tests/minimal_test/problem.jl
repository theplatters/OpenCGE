# ── Problem function (mobile labor CES with numeraire) ──

function mobile_labor_problem(X, params)
    (; Ω, consumption_share, factor_share, supply_shock, demand_shock, labor_bar, θ, ϵ, σ, η) = params
    N = length(factor_share)

    p = max.(X[1:N], 1e-10)
    y = max.(X[N+1:2N], 1e-10)
    w = max(X[2N+1], 1e-10)

    out = zeros(eltype(X), 2N + 1)

    # Intermediate goods price index
    intermediate_price = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

    # CPI
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Sectoral labor demand
    L_i = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ

    # Final demand
    total_income = w * sum(L_i)
    final_demand = (total_income * p .^ (-σ) .* demand_shock .* consumption_share) ./ cpi .^ (-σ)

    # Intermediary demand
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # Cost function
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    # Eq 1: Price equations for sectors 2..N (drop sector 1 → numeraire)
    out[1:N-1] = p[2:N] .- cost[2:N]

    # Eq 2: Market clearing (N equations)
    out[N:2N-1] = y .- intermediary_demand .- final_demand

    # Eq 3: Labor market clearing
    out[2N] = sum(L_i) - labor_bar * w^η

    # Eq 4: Numeraire constraint — CPI = 1
    out[2N+1] = cpi - 1.0

    return out
end