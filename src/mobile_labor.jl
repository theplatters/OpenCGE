# ═══════════════════════════════════════════════════════════════════════════════
# Mobile Labor with Continuous Labor Supply Elasticity η
# ═══════════════════════════════════════════════════════════════════════════════
#
# This file implements the key model extension: replacing sector-specific
# (immobile) labor with mobile labor that responds to an economy-wide wage
# through a continuous elasticity parameter η ∈ [0, ∞).
#
#   η = 0     → perfectly inelastic labor (fixed total supply, like Leontief closure)
#   η → ∞     → perfectly elastic labor (full slack, unlimited workers)
#   0 < η < ∞ → intermediate: L = L̄ · w^η
#
# The key change vs. the base CES model:
#   - The unknown vector gains a scalar wage w: X = [p(1:N); y(1:N); w]
#   - Labor allocation is endogenous: L_i = (∂Y_i/∂L_i = w) → solved from FOC
#   - Labor market clearing: Σ L_i(p,y,w) = L̄ · w^η  (the new equation)
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    MobileLaborCESElasticities

Elasticity parameters for the mobile-labor CES model.

- `θ` : elasticity of substitution between intermediate goods
- `ϵ`  : elasticity of substitution between labor and intermediate composite
- `σ`  : elasticity of substitution in consumption
- `η`  : labor supply elasticity (0 = inelastic, ∞ = perfectly elastic)
"""
struct MobileLaborCESElasticities <: AbstractElasticities
    θ::Float64
    ϵ::Float64
    σ::Float64
    η::Float64
end

"""
    MobileLaborCES

Model type for CES with mobile labor and endogenous labor supply elasticity.

- `elasticities` : a `MobileLaborCESElasticities` struct
- `labor_bar`    : baseline total labor supply L̄ (defaults to Σ labor_share)
"""
struct MobileLaborCES <: ModelType
    elasticities::MobileLaborCESElasticities
    labor_bar::Float64
end

function MobileLaborCES(elasticities::MobileLaborCESElasticities, data::Data)
    labor_bar = sum(data.labor_share)
    MobileLaborCES(elasticities, labor_bar)
end

# ─────────────────────────────────────────────────────────────────────────────────
# Sectoral labor demand from the wage (marginal product condition)
# ─────────────────────────────────────────────────────────────────────────────────
# In the base CES model, the wage equation (from FOC of labor) is:
#   w_i = p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) · L_i^(-1/ε)
#
# With mobile labor, w is economy-wide (scalar). Inverting for L_i:
#   L_i = [ p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) / w ]^ε
#
# This gives sectoral labor demand as a function of (p, y, w).
# ─────────────────────────────────────────────────────────────────────────────────

"""
    sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})

Compute sectoral labor demand L_i given prices, quantities, and the economy-wide wage.
Derived by inverting the marginal-product-of-labor condition.
"""
function sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; supply_shock) = shocks
    (; factor_share) = data

    # L_i = [ p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) / w ]^ε
    return (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ
end

"""
    economy_wide_wage(p, y, labor, model::Model{MobileLaborCES})

Compute the sector-implied wage from the base CES FOC (used for initialization).
When labor is mobile, all sectors should imply the same wage w.
"""
function economy_wide_wage(p, y, labor, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; supply_shock) = shocks
    (; factor_share) = data

    # w_i = p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) · L_i^(-1/ε)
    w_vec = p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* (labor .^ (-1 / ϵ))
    # Return the GDP-weighted average as a single wage
    return sum(w_vec .* labor) / sum(labor)
end

# ═══════════════════════════════════════════════════════════════════════════════
# The equilibrium problem (2N+1 equations, 2N+1 unknowns)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    problem(out, X, model::Model{MobileLaborCES})

The equilibrium system for the mobile-labor CES model.

Unknowns: X = [p(1:N); y(1:N); w]  — 2N+1 elements
Equations (2N+1):
  1. Price equations        (N-1):  p_i = cost_i(p, w)   for i=2..N
  2. Market clearing        (N):     y_i = intermediary_demand_i + final_demand_i
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄ · w^η
  4. Numeraire              (1):     CPI = 1  (Σ β_i · p_i^(1-σ))^(1/(1-σ)) = 1)

Note: The first price equation (i=1) is replaced by the numeraire constraint.
This breaks the price-level indeterminacy inherent in CRTS models.

Economic note: Under CRTS + mobile labor, the zero-profit condition pins down
the wage w and relative prices independent of η. The labor supply elasticity η
determines the scale of the economy (total employment, GDP) but not relative
prices. This is a standard result in CGE modeling.
"""
function problem(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)

    # Unpack unknowns
    p = max.(X[1:N], 0)
    y = max.(X[N+1:2N], 0)
    w = max(X[2N+1], 1e-10)  # scalar wage, keep positive

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # ── Intermediate goods price index ──
    intermediate_price = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

    # ── CPI (consumption price index) ──
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # ── Sectoral labor demand (from inverted FOC) ──
    L_i = sectoral_labor_demand(p, y, w, model)

    # ── Final demand ──
    total_income = w * sum(L_i)
    final_demand = (total_income * p .^ (-σ) .* demand_shock .* consumption_share) ./ cpi .^ (-σ)

    # ── Intermediary demand ──
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # ── Cost function ──
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    # ── Equation 1: Price equations for sectors 2..N (N-1 equations) ──
    out[1:N-1] .= p[2:N] .- cost[2:N]

    # ── Equation 2: Market clearing (N equations) ──
    out[N:2N-1] .= y .- intermediary_demand .- final_demand

    # ── Equation 3: Labor market clearing (1 equation) ──
    out[2N] = sum(L_i) - options.labor_bar * w^η

    # ── Equation 4: Numeraire constraint — CPI = 1 ──
    out[2N+1] = cpi - 1.0

    nothing
end

# ═══════════════════════════════════════════════════════════════════════════════
# Solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    solve(model::Model{MobileLaborCES}; init)

Solve the mobile-labor CES model. Returns a `Solution` with the equilibrium
prices, quantities, wage (as a vector of the same wage in all sectors),
consumption, and GDP measures.

The initial guess defaults to baseline prices=1, quantities=λ, and wage=1.
"""
function solve(model::Model{MobileLaborCES};
    init = nothing
)
    (; data, options, shocks) = model
    N = length(data.factor_share)

    if init === nothing
        # Default initialization: p=1, y=λ, w=1
        init = [ones(N); data.λ; 1.0]
    end

    ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
    x = NonlinearSolve.solve(ProbN, reltol=1e-8, abstol=1e-8).u

    p = x[1:N]
    q = x[N+1:2N]
    w = x[2N+1]

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all equal to w — mobile labor)
    wages = fill(w, N)

    # Consumption
    consumption_share_adj = shocks.demand_shock .* data.consumption_share
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    total_income = w * sum(L_i)
    consumption = total_income .* consumption_share_adj .* (p / numeraire) .^ (-σ)

    # GDP measures
    laspeyres_index = sum(consumption) / sum(data.consumption_share)
    nominal_gdp = (w * sum(L_i)) / numeraire

    return Solution(p, q, wages, consumption, numeraire, laspeyres_index, nominal_gdp, model)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility: construct a model from data + elasticities
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=nothing)

Convenience constructor for a MobileLaborCES model.
"""
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Float64, Nothing}=nothing)
    el = MobileLaborCESElasticities(θ, ϵ, σ, η)
    lb = labor_bar === nothing ? sum(data.labor_share) : labor_bar
    model = Model(data, shocks, MobileLaborCES(el, lb))
    return model
end
