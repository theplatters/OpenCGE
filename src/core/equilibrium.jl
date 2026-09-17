# src/core/equilibrium.jl — Solution + mobile-labor equilibrium
# (concatenated verbatim: src/solution.jl + src/mobile_labor.jl, in that order;
# _closure_symbol moved to src/closures/labor/types.jl)

struct Solution
	prices::Vector{Float64}
	prices_raw::Vector{Float64}
	quantities::Vector{Float64}
	wages::Vector{Float64}
	wages_raw::Vector{Float64}
	consumption::Vector{Float64}
	numeraire::Float64
	real_gdp::Float64
	nominal_gdp::Float64
	model::Model
end

function Solution(prices_raw, quantities, wages, consumption, numeraire, real_gdp, nominal_gdp, model)
	return Solution(
		prices_raw ./ numeraire,
		prices_raw,
		quantities,
		wages ./ numeraire,
		wages,
		consumption,
		numeraire,
		real_gdp,
		nominal_gdp,
		model)
end

function real_gdp(sol::Solution)::Float64
	return sol.real_gdp
end

function nominal_gdp(sol::Solution)::Float64
	return sol.nominal_gdp
end

"""Return exact equilibrium residuals at the raw values stored in `sol`."""
function equilibrium_residuals(sol::Solution)
	_equilibrium_residuals(sol.model, sol)
end

_equilibrium_residuals(::Model, ::Solution) =
	throw(ArgumentError("equilibrium_residuals is not supported for this model type"))

"""
	tornqvist_quantity_index(prices, quantities, base_prices, base_quantities)

Compute the discrete Törnqvist approximation to the Divisia quantity index.
The logarithmic growth rate is the expenditure-share-weighted change in final
quantities,

`log(Q₁/Q₀) = sumᵢ 0.5 * (sᵢ₀ + sᵢ₁) * log(qᵢ₁/qᵢ₀)`.

Goods with zero expenditure shares in both equilibria are omitted. A good with
a positive expenditure share in either equilibrium must have a strictly
positive quantity in both equilibria.
"""
function tornqvist_quantity_index(
	prices::AbstractVector{<:Real},
	quantities::AbstractVector{<:Real},
	base_prices::AbstractVector{<:Real},
	base_quantities::AbstractVector{<:Real},
)::Float64
	lengths = length.((prices, quantities, base_prices, base_quantities))
	all(==(first(lengths)), lengths) ||
		throw(DimensionMismatch("prices and quantities must have matching lengths"))
	any(x -> x < 0, prices) && throw(DomainError(prices, "prices must be nonnegative"))
	any(x -> x < 0, base_prices) && throw(DomainError(base_prices, "base prices must be nonnegative"))
	any(x -> x < 0, quantities) && throw(DomainError(quantities, "quantities must be nonnegative"))
	any(x -> x < 0, base_quantities) && throw(DomainError(base_quantities, "base quantities must be nonnegative"))

	current_expenditure = prices .* quantities
	base_expenditure = base_prices .* base_quantities
	current_total = sum(current_expenditure)
	base_total = sum(base_expenditure)
	current_total > 0 || throw(DomainError(current_total, "current expenditure must be positive"))
	base_total > 0 || throw(DomainError(base_total, "base expenditure must be positive"))

	current_shares = current_expenditure ./ current_total
	base_shares = base_expenditure ./ base_total
	active = (current_shares .+ base_shares) .> 0
	any(quantities[active] .<= 0) &&
		throw(DomainError(quantities, "positive-share goods require positive quantities"))
	any(base_quantities[active] .<= 0) &&
		throw(DomainError(base_quantities, "positive-share goods require positive base quantities"))

	log_growth = sum(
		0.5 .* (base_shares[active] .+ current_shares[active]) .*
		log.(quantities[active] ./ base_quantities[active]),
	)
	return exp(log_growth)
end

function wages(sol::Solution)::Vector{Float64}
	return sol.wages
end

function consumption(sol::Solution)::Vector{Float64}
	sol.consumption
end

function Base.getindex(sol::Solution, ::Colon, sector::String)
	index = findfirst(==(sector), sol.model.data.io.Sektoren)
	return Dict(:prices => sol.prices[index], :quantities => sol.quantities[index], :wages => sol.wages[index], :consumption => sol.consumption[index])
end

struct SectorData
	name::String
	price::Float64
	quantity::Float64
	wage::Float64
	consumption::Float64
end

function eachsector(sol::Solution)
	sectors = sol.model.data.io.Sektoren
	return (SectorData(sectors[i], sol.prices[i], sol.quantities[i], sol.wages[i], sol.consumption[i]) for i in 1:length(sol.prices))
end

function Base.getindex(sol::Solution, sector_index::Int)
	sectors = sol.model.data.io.Sektoren
	return SectorData(
		sectors[sector_index],
		sol.prices[sector_index],
		sol.quantities[sector_index],
		sol.wages[sector_index],
		sol.consumption[sector_index],
	)
end
function Base.getindex(sol::Solution, sector_indices::UnitRange{Int})
	sectors = sol.model.data.io.Sektoren
	return [SectorData(
		sectors[i],
		sol.prices[i],
		sol.quantities[i],
		sol.wages[i],
		sol.consumption[i],
	) for i in sector_indices]
end
function Base.getindex(sol::Solution, indices::Vector{Int})
	sectors = sol.model.data.io.Sektoren
	return [SectorData(
		sectors[i],
		sol.prices[i],
		sol.quantities[i],
		sol.wages[i],
		sol.consumption[i],
	) for i in indices]
end

function multiplier(sol::Solution)::Float64
	(; data, shocks) = sol.model
	simple_effect = 1 + sum(shocks.demand_shock_raw) ./ sum(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])
	(1 .- sol.real_gdp) ./ (1 .- simple_effect)
end

# --- src/mobile_labor.jl (verbatim modulo the moved _closure_symbol block) ---

# ═══════════════════════════════════════════════════════════════════════════════
# Mobile Labor with Geometric Intersectoral Reallocation η
# ═══════════════════════════════════════════════════════════════════════════════
#
# This file implements the key model extension: replacing sector-specific
# (immobile) labor with geometric intersectoral reallocation parameter η.
#
#   η = 0     → immobile baseline allocation
#   η = 1     → fully cost-minimizing allocation
#   η outside [0, 1] → extrapolation of the same geometric rule
#
# The key change vs. the base CES model:
#   - The unknown vector gains a scalar wage w: X = [p(1:N); y(1:N); w]
#   - Labor allocation is endogenous: L_i = (∂Y_i/∂L_i = w) → solved from FOC
#   - Labor market clearing: total employment equals the fixed labor bar
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    MobileLaborCESElasticities

Elasticity parameters for the mobile-labor CES model.

- `θ` : elasticity of substitution between intermediate goods
- `ϵ`  : elasticity of substitution between labor and intermediate composite
- `σ`  : elasticity of substitution in consumption
- `η`  : intersectoral reallocation parameter; it is not a labor-supply elasticity
"""
struct MobileLaborCESElasticities <: AbstractElasticities
    θ::Float64
    ϵ::Float64
    σ::Float64
    η::Float64
end

"""
    MobileLaborCES

Model type for CES with geometric intersectoral labor reallocation and an explicit wage-regime closure.

- `elasticities` : a `MobileLaborCESElasticities` struct
- `labor_bar`    : fixed total employment capacity L̄ (defaults to Σ labor_share)
"""
struct MobileLaborCES <: ModelType
    elasticities::MobileLaborCESElasticities
    labor_bar::Float64
    closure::Symbol   # :mobile = flexible wage; :fixed = sticky wage
    function MobileLaborCES(elasticities::MobileLaborCESElasticities, labor_bar::Float64, closure::Symbol)
        closure in (:mobile, :fixed) || throw(ArgumentError("unsupported MobileLaborCES closure $closure; use :mobile or :fixed"))
        new(elasticities, labor_bar, closure)
    end
end

# _closure_symbol lives in src/closures/labor/types.jl (moved; supports :beta as well).

MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::AbstractLaborClosure) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real; closure=:mobile) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::Symbol) = MobileLaborCES(e, Float64(lb), closure)

function MobileLaborCES(elasticities::MobileLaborCESElasticities, data::Data)
    labor_bar = sum(data.labor_share)
    MobileLaborCES(elasticities, labor_bar, :mobile)
end

labor_closure(options::MobileLaborCES) = options.closure == :mobile ? FlexibleWageClosure() : FixedWageClosure()

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

const _ETA_MAX = 50.0
const _ETA_SCALE_INDETERMINACY_TOL = 1.1e-6
const _ALLOCATION_LOG_LIMIT = 50.0
const _WEDGE_MIN = eps(Float64)

function _checked_eta(η)
    isfinite(η) || throw(ArgumentError("η must be finite"))
    abs(η) <= _ETA_MAX || throw(DomainError(η, "η must have magnitude at most 50"))
    Float64(η)
end

_positive_floor(x) = max.(x, eps(Float64))
_positive_floor(x::Real) = max(x, eps(Float64))

"""Cost-minimizing labor demand, evaluated in log space for stable η extrapolation."""
function _cost_minimizing_labor(p, y, w, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; factor_share) = data
    p, y, w, A, α = _positive_floor.((p, y, w, shocks.supply_shock, factor_share))
    log_demand = ϵ .* (log.(p) .+ ((ϵ - 1) / ϵ) .* log.(A) .+ (1 / ϵ) .* log.(α) .+
        (1 / ϵ) .* log.(y) .- log(w))
    exp.(clamp.(log_demand, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""Interpolate fixed and cost-minimizing allocations geometrically in log space."""
function _interpolated_labor(p, y, w, model::Model{MobileLaborCES})
    η = _checked_eta(model.options.elasticities.η)
    fixed = _positive_floor(model.data.labor_share)
    optimum = _positive_floor(_cost_minimizing_labor(p, y, w, model))
    log_labor = (1 - η) .* log.(fixed) .+ η .* log.(optimum)
    exp.(clamp.(log_labor, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""Return the allocative-efficiency factor (one at the optimum, never above one).

The log allocation ratio is limited to ±50 to keep trial solver iterates finite.
"""
function _allocation_efficiency_wedge(labor, optimum, factor_share, ϵ)
    ϵ > 0 || throw(DomainError(ϵ, "ϵ must be positive for the allocation wedge"))
    ratio = clamp.(log.(_positive_floor(labor) ./ _positive_floor(optimum)),
                   -_ALLOCATION_LOG_LIMIT, _ALLOCATION_LOG_LIMIT)
    curvature = abs(1 - ϵ) / ϵ
    penalty = 0.5 .* _positive_floor(factor_share) .* max.(1 .- factor_share, 0) .* curvature .* ratio.^2
    # Keep the wedge strictly positive even when the clamped log ratio makes
    # the exponential underflow. `max.` remains elementwise AD-compatible.
    max.(exp.(-penalty), _WEDGE_MIN)
end

function _allocation_efficiency_wedge(p, y, w, labor, model::Model{MobileLaborCES})
    (; factor_share) = model.data
    ϵ = model.options.elasticities.ϵ
    optimum = _cost_minimizing_labor(p, y, w, model)
    _allocation_efficiency_wedge(labor, optimum, factor_share, ϵ)
end

"""
    sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})

Compute sectoral labor demand L_i given prices, quantities, and the economy-wide wage.
Derived by inverting the marginal-product-of-labor condition.
"""
function sectoral_labor_demand(p, y, w, model::Model{MobileLaborCES})
    _interpolated_labor(p, y, w, model)
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
  1. Zero-profit equations  (N):     p_i = cost_i(p, w)   for all i=1..N
  2. Market clearing        (N-1):   y_i = intermediary_demand_i + final_demand_i
                                      for i=1..N-1 (the last sector's equation is
                                      dropped under Walras' law)
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄
  4. Numeraire              (1):     CPI = 1  (Σ β_i · p_i^(1-σ))^(1/(1-σ)) = 1)

Note: The numeraire constraint (CPI = 1) pins the price level, breaking the
price-level indeterminacy inherent in CRTS models. All N zero-profit
conditions are enforced because the numeraire already replaces the price-level
degree of freedom.

Economic note: η changes the geometric allocation between baseline and
cost-minimizing sectoral labor demand. It is not a labor-supply elasticity.
"""
function problem(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)

    # Unpack unknowns
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)  # scalar wage, keep positive

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # ── Intermediate goods price index ──
    intermediate_price = (Ω_raw * p .^ (1 - θ)) .^ (1 / (1 - θ))

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
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # ── B&F (2019) labor-reallocation wedge ──
    # L_opt is the cost-minimizing (mobile) labor demand; L_base is the
    # baseline (immobile, fixed) allocation. When labor cannot reallocate toward
    # its shocked optimum (eta < 1), there is a second-order allocative
    # inefficiency that raises effective cost / lowers TFP. The wedge vanishes
    # for eta=1 (mobile) and at baseline (L_opt == L_base).
    alloc_wedge = _allocation_efficiency_wedge(p, y, w, L_i, model)

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

    # This is exclusively the flexible-wage system; fixed wages are routed to
    # `problem_fixed` by `solve`.
    out[2N] = sum(L_i) - options.labor_bar

    # ── Equation 4: Numeraire constraint — CPI = 1 ──
    out[2N+1] = cpi - 1.0

    nothing
end

"""Return the exact residual vector for either mobile-labor closure."""
function equilibrium_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    fixed = labor_closure(model.options) isa FixedWageClosure
    expected = fixed ? 2N : 2N + 1
    length(X) == expected || throw(DimensionMismatch("closure expects a $expected element vector"))
    out = zeros(Float64, expected)
    if fixed
        problem_fixed(out, collect(X), model)
    else
        problem(out, collect(X), model)
    end
    out
end

function _equilibrium_residuals(model::Model{MobileLaborCES}, sol::Solution)
    X = labor_closure(model) isa FixedWageClosure ?
        [sol.prices_raw; sol.quantities] :
        [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
    equilibrium_residuals(model, X)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Solver
# ═══════════════════════════════════════════════════════════════════════════════

"""
    problem_fixed(out, X, model)

Specialized equilibrium system for the :fixed (sticky-wage) closure.
w = 1.0 is hard-coded; unknowns are p(1:N) and y(1:N) only (2N).
Equations:
  1. N zero-profit conditions:  p_i = cost_i(p, w=1.0)
  2. N-1 market-clearing:       y_i = intermed_i + final_i  (drop last sector)
  3. 1 numeraire:               CPI = 1  (replaces the redundant market equation)
Total: 2N equations, 2N unknowns.
Employment `L_i` is computed post-solve and is not constrained to `labor_bar`.
The sticky-wage counterfactual holds the production cost at its direct CES
form. The mobile-only reallocation-efficiency correction is not added because
this fixed closure has no labor-market/accounting equation to support it.
"""
function problem_fixed(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)
    w = 1.0  # sticky wage

    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # Intermediate goods price index
    intermediate_price = (Ω_raw * p .^ (1 - θ)) .^ (1 / (1 - θ))

    # CPI
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Sectoral labor demand at w=1.0
    L_i = sectoral_labor_demand(p, y, w, model)

    # Final demand (budget-consistent)
    total_income = w * sum(L_i)
    agg = sum(consumption_share .* demand_shock .* p .^ (1 - σ))
    final_demand = (consumption_share .* demand_shock .* total_income .* p .^ (-σ)) ./ agg

    # Autonomous & investment final demand
    cons_base = sum(labor_share)
    A = shocks.autonomous_demand .* consumption_share .* cons_base
    G = shocks.investment_shock .* consumption_share .* cons_base
    total_final_demand = final_demand .+ A .+ G

    # Intermediary demand
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # Direct CES cost at the sticky wage; the reallocation correction is
    # supported only by the mobile closure's labor-market accounting.
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    # 1. Zero-profit (all N sectors)
    out[1:N] .= p .- cost

    # 2. Market clearing (N-1 equations, drop last)
    out[N+1:2N-1] .= y[1:N-1] .- intermediary_demand[1:N-1] .- total_final_demand[1:N-1]

    # 3. Numeraire: CPI = 1 (2N-th equation)
    out[2N] = cpi - 1.0

    nothing
end


"""
    _solve_fixed(model; init)

Solve the sticky-wage (w=1.0) system. Returns a `Solution` with equilibrium
prices, quantities, wage=1.0, consumption, and Tornqvist real GDP computed
from consumption (B&F metric). Employment (`sum(L_i)`) is a post-solve outcome.
"""
function _solve_fixed(model::Model{MobileLaborCES}; init=nothing)
    (; data, options, shocks) = model
    N = length(data.factor_share)

    η = options.elasticities.η
    isapprox(η, 1.0; rtol=0, atol=_ETA_SCALE_INDETERMINACY_TOL) &&
        all(iszero, shocks.autonomous_demand) && all(iszero, shocks.investment_shock) &&
        throw(ArgumentError("fixed-wage η=1 has a homogeneous, scale-indeterminate equilibrium; add autonomous or investment demand as an additive-demand anchor, or use another η"))

    if init === nothing
        init = [ones(N); data.λ]
    elseif length(init) == 2N + 1
        # Drop the wage component if a :mobile init was provided
        init = init[1:2N]
    end

    # Avoid asking the nonlinear solver to differentiate an already exact
    # baseline, which some solver versions report as stalled.
    x = if maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
        Float64.(init)
    else
        ProbN = NonlinearSolve.NonlinearProblem(problem_fixed, init, model)
        res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=1000)
        string(res.retcode) == "Success" ||
            error("MobileLaborCES._solve_fixed did not converge: retcode = $(res.retcode)")
        res.u
    end

    p = x[1:N]
    q = x[N+1:2N]
    w = 1.0  # sticky wage

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all 1.0 — sticky)
    wages = fill(w, N)

    # Consumption (budget-consistent)
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    total_income = w * sum(L_i)
    agg = sum(data.consumption_share .* shocks.demand_shock .* p .^ (1 - σ))
    consumption = (data.consumption_share .* shocks.demand_shock .* total_income .* p .^ (-σ)) ./ agg

    # Real GDP: consumption Tornqvist (B&F metric)
    base_income = sum(data.labor_share)
    base_consumption = data.consumption_share .* base_income
    real_gdp_index = tornqvist_quantity_index(
        p,
        consumption,
        ones(N),
        base_consumption,
    )
    nominal_gdp = w * sum(L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end


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

    # ── Sticky-wage path (:fixed closure) ──
    # w = 1.0 is hard-coded; solve a 2N system (prices + quantities only).
    # Employment is computed post-solve.
    if labor_closure(options) isa FixedWageClosure
        return _solve_fixed(model; init=init)
    end

    # ── Standard :mobile path (2N+1 with wage as unknown) ──
    if init === nothing
        # Default initialization: p=1, y=λ, w=1
        init = [ones(N); data.λ; 1.0]
    end

    # As in the fixed-wage path, avoid asking the nonlinear solver to
    # differentiate an already exact baseline (some versions report this as
    # stalled). Non-exact guesses still retain retcode checks.
    x = if maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
        Float64.(init)
    else
        ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
        res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=1000)
        # `retcode` may be a symbol or string depending on NonlinearSolve.
        string(res.retcode) == "Success" ||
            error("MobileLaborCES.solve did not converge: retcode = $(res.retcode)")
        res.u
    end

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
    nominal_gdp = w * sum(L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility: construct a model from data + elasticities
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=nothing)

Convenience constructor for a MobileLaborCES model.
`closure=:mobile` uses a flexible, market-clearing wage; `closure=:fixed` uses a
sticky wage of one and unconstrained employment demand.
"""
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Real, Nothing}=nothing, closure=:mobile)
    el = MobileLaborCESElasticities(θ, ϵ, σ, η)
    closure_symbol = _closure_symbol(closure)
    labor_bar !== nothing && closure_symbol == :fixed && throw(ArgumentError(
        "fixed closure treats employment as an outcome; labor_bar is not used and must not be supplied"))
    lb = labor_bar === nothing ? sum(data.labor_share) : Float64(labor_bar)
    model = Model(data, shocks, MobileLaborCES(el, lb, closure_symbol))
    return model
end

mobile_labor_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real;
    labor_bar=nothing, closure=:mobile) = mobile_labor_model(data, shocks, Float64(θ), Float64(ϵ), Float64(σ), Float64(η);
    labor_bar=labor_bar === nothing ? nothing : Float64(labor_bar),
    closure=_closure_symbol(closure))
