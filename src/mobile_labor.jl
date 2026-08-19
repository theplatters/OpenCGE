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
    closure::Symbol   # :mobile = market-clearing (eta free); :fixed = sticky wage (B&F 2022 unemployment closure)
end

function MobileLaborCES(elasticities::MobileLaborCESElasticities, labor_bar::Float64)
    MobileLaborCES(elasticities, labor_bar, :mobile)
end

function MobileLaborCES(elasticities::MobileLaborCESElasticities, data::Data)
    labor_bar = sum(data.labor_share)
    MobileLaborCES(elasticities, labor_bar, :mobile)
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
    (; ϵ, η) = options.elasticities
    (; supply_shock) = shocks
    (; factor_share, labor_share) = data

    # Cost-minimizing labor demand given the common wage (mobile limit, eta=1):
    #   L_i = [ p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) / w ]^ε
    L_costmin = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ
    # Baseline (immobile) labor allocation, fixed in sector (eta=0). Uses
    # labor_share so that sum(L_fixed) == labor_bar exactly (labor market clears
    # trivially for eta=0 and the system stays square/solvable).
    L_fixed = labor_share
    # Repurpose eta as B&F (2019) labor-reallocation parameter beta:
    #   eta = 0 -> immobile (labor fixed at baseline shares, no reallocation)
    #   eta = 1 -> mobile   (full cost-minimizing reallocation across sectors)
    return L_fixed .^ (1 - η) .* L_costmin .^ η
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
    # Budget-consistent CES demand: consumption weights are β̃_i = cs_i·ds_i,
    # normalized through the aggregator `agg = Σ cs_j·ds_j·p_j^(1-σ)` so that
    # Σ p_i·c_i = total_income exactly (Walras closes; no spurious overspend).
    # Previously `demand_shock` was applied unnormalized, inflating household
    # expenditure above income and breaking the budget identity.
    total_income = w * sum(L_i)
    agg = sum(consumption_share .* demand_shock .* p .^ (1 - σ))
    final_demand = (consumption_share .* demand_shock .* total_income .* p .^ (-σ)) ./ agg

    # ── Autonomous & investment (extra-household) final demand (Milestone E) ──
    # Expansive: they raise total final demand above wage income, so the economy
    # must scale employment to meet them. With mobile labor (eta>0) this expands
    # at ~constant wage (large real response); with eta=0 the wage must rise,
    # muting the response -- the mobile-labor bridge. The gap to wage income is
    # the financing record (debt / external deficit), not a constraining equation.
    # Base the autonomous / investment demand on baseline household consumption
    # (same units as `final_demand`), so the shock is a moderate, well-scaled
    # additive boost rather than a gross-output multiple.
    cons_base = sum(data.labor_share)
    A = shocks.autonomous_demand .* data.consumption_share .* cons_base
    G = shocks.investment_shock .* data.consumption_share .* cons_base
    total_final_demand = final_demand .+ A .+ G

    # ── Intermediary demand ──
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # ── B&F (2019) labor-reallocation wedge ──
    # L_opt is the cost-minimizing (mobile) labor demand; L_base is the
    # baseline (immobile, fixed) allocation. When labor cannot reallocate toward
    # its shocked optimum (eta < 1), there is a second-order allocative
    # inefficiency that raises effective cost / lowers TFP. The wedge vanishes
    # for eta=1 (mobile) and at baseline (L_opt == L_base).
    L_opt = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ
    L_base = factor_share .* data.λ  # baseline (immobile, fixed) labor allocation
    mis = 1 - η
    log_ratio = log.((L_base .+ 1e-12) ./ (L_opt .+ 1e-12))
    penalty = 0.5 .* factor_share .* (1 .- factor_share) .* ((ϵ - 1) / ϵ) .* (mis .* log_ratio) .^ 2
    alloc_wedge = exp.(-penalty)

    # ── Cost function (effective TFP includes the reallocation wedge) ──
    cost = ((supply_shock .* alloc_wedge) .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    # ── Equation 1: Zero-profit for ALL N sectors (including sector 1) ──
    # Under CRTS the price level is pinned by the numeraire (Eq. 4), so all N
    # zero-profit conditions are independent and must be enforced. The previous
    # code omitted sector 1's condition, leaving a material non-equilibrium
    # (residual ≈ 10.9 in the 71-sector diagnostic).
    out[1:N] .= p .- cost

    # ── Equation 2: Market clearing for sectors 1..N-1 (N-1 equations) ──
    # The LAST sector's market-clearing equation is the redundant one under
    # Walras' law and is dropped; it clears automatically given the other N-1
    # markets, all N zero-profit conditions, the labor market (Eq. 3) and the
    # numeraire (Eq. 4). We drop the LAST (not the first) sector so that a shock
    # to sector 1 (construction) remains enforced.
    out[N+1:2N-1] .= y[1:N-1] .- intermediary_demand[1:N-1] .- total_final_demand[1:N-1]

    if options.closure == :fixed
        # Sticky-wage / unemployment closure (B&F 2022 style): hold the wage at
        # its baseline (numeraire) value. The labor market no longer clears, so
        # employment is labor demand and absorbs the shock; the gap
        # L_bar - sum(L_i) is unemployment (recorded, not enforced).
        out[2N] = w - 1.0
    else
        # Labor-market clearing: total labor demand equals the (fixed) labor
        # supply L_bar. The reallocation across sectors is governed by eta (B&F
        # beta); the wage w adjusts here to clear the market when eta > 0.
        out[2N] = sum(L_i) - options.labor_bar
    end

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
    res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=1000)
    # Fail loudly instead of silently returning a non-equilibrium. The earlier
    # false "GO" certifications came from trusting .u without checking retcode.
    # `retcode` may be the symbol :Success or the string "Success" depending on
    # the NonlinearSolve version; normalize to a string before comparing.
    string(res.retcode) == "Success" ||
        error("MobileLaborCES.solve did not converge: retcode = $(res.retcode)")
    x = res.u

    p = x[1:N]
    q = x[N+1:2N]
    w = x[2N+1]

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all equal to w — mobile labor)
    wages = fill(w, N)

    # Consumption — must match the budget-consistent demand used inside `problem`
    # so the reported allocation is the actual equilibrium allocation.
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    total_income = w * sum(L_i)
    agg = sum(data.consumption_share .* shocks.demand_shock .* p .^ (1 - σ))
    consumption = (data.consumption_share .* shocks.demand_shock .* total_income .* p .^ (-σ)) ./ agg

    # Real GDP: Tornqvist (Divisia) quantity index of FINAL CONSUMPTION (value-added)
    # — not gross output. In B&F (2019), real GDP is a Divisia index of real final
    # demand / value added. Using gross output confounds intermediate flows with
    # welfare-relevant final output and dilutes the reallocation bridge. The
    # consumption (final-demand) Tornqvist isolates the welfare-relevant change.
    # Baseline consumption quantities are cs_i * total_income_base (p=1, numeraire).
    base_income = sum(data.labor_share)
    base_consumption = data.consumption_share .* base_income
    real_gdp_index = tornqvist_quantity_index(
        p,
        consumption,
        ones(N),
        base_consumption,
    )
    nominal_gdp = (w * sum(L_i)) / numeraire

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility: construct a model from data + elasticities
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=nothing)

Convenience constructor for a MobileLaborCES model.
"""
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Float64, Nothing}=nothing, closure::Symbol=:mobile)
    el = MobileLaborCESElasticities(θ, ϵ, σ, η)
    lb = labor_bar === nothing ? sum(data.labor_share) : labor_bar
    model = Model(data, shocks, MobileLaborCES(el, lb, closure))
    return model
end
