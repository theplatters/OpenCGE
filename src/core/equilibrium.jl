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

# --- src/mobile_labor.jl (ported from cbase2/src/core/mobile_labor.jl with ADR-0005 compatibility (manna retained, legacy error surfaces)) ---

# Backport of the recorded kernel DIFFs (ADR-0005): eta_s/:beta hook,
# financing hooks, _intermediate_price/_ces_unit_cost guards, residual-based
# solver gates with LM polish, all-N fixed formulation, household_baseline
# Törnqvist base, financing/eta_s kwargs. Compatibility: legacy manna A/G
# retained; legacy error substrings kept.
# Port of the cbase2 review fixes (ADR-0010, 2026-09-17): the MOBILE system
# keeps N-1 clearing equations plus the CPI = 1 numeraire and exposes the
# omitted N-th market as the residual external account
# (`market_clearing_residuals`); the allocation wedge is retired with η ∈ {0,1}
# (the all-N `644ba37` form, later reverted on `revisefinal` by `0f33ad6`,
# over-determines the open economy and has no exact root for additive-demand
# cases). See `problem`.

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
    # BETA closure (Stage 1.2): elasticity of TOTAL labour supply along the
    # real wage; 0.0 = vertical supply (ALPHA). Used only by the :beta closure.
    eta_s::Float64
end
# Parent-compatible 4-arg constructor (no elastic labour supply).
MobileLaborCESElasticities(θ::Real, ϵ::Real, σ::Real, η::Real) =
    MobileLaborCESElasticities(Float64(θ), Float64(ϵ), Float64(σ), Float64(η), 0.0)

"""
    MobileLaborCES

Model type for CES with geometric intersectoral labor reallocation and an explicit wage-regime closure.

- `elasticities` : a `MobileLaborCESElasticities` struct
- `labor_bar`    : fixed total employment capacity L̄ (defaults to Σ labor_share)
"""
struct MobileLaborCES <: ModelType
    elasticities::MobileLaborCESElasticities
    labor_bar::Float64
    closure::Symbol   # :mobile = flexible wage (ALPHA); :fixed = sticky wage (GAMMA);
                      # :beta = flexible wage with elastic total labour supply (BETA)
    function MobileLaborCES(elasticities::MobileLaborCESElasticities, labor_bar::Float64, closure::Symbol)
        closure in (:mobile, :fixed, :beta) || throw(ArgumentError("unsupported MobileLaborCES closure $closure; use :mobile, :fixed, or :beta"))
        new(elasticities, labor_bar, closure)
    end
end

MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::AbstractLaborClosure) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real; closure=:mobile) =
    MobileLaborCES(e, Float64(lb), _closure_symbol(closure))
MobileLaborCES(e::MobileLaborCESElasticities, lb::Real, closure::Symbol) = MobileLaborCES(e, Float64(lb), closure)

function MobileLaborCES(elasticities::MobileLaborCESElasticities, data::Data)
    labor_bar = sum(data.labor_share)
    MobileLaborCES(elasticities, labor_bar, :mobile)
end

labor_closure(options::MobileLaborCES) =
    options.closure == :mobile ? FlexibleWageClosure() :
    options.closure == :beta   ? ElasticLaborClosure(options.elasticities.eta_s) :
                                 FixedWageClosure()

# ── Labor-market equation hook ──
# The total-labor-market residual of the flexible-wage systems (problem(), 2N+1).
# ALPHA (FlexibleWageClosure): vertical supply at the bar L̄ (the CPI argument is
#       unused).
# BETA  (ElasticLaborClosure, defined in src/closures/labor/types.jl): elastic
#       supply on the REAL wage, L^s = L̄ · [(w/P)/(w0/P0)]^{η_s}, deflated by the
#       passed CPI (ADR-0014; DE-0004 requires the real wage).
labor_market_residual(::FlexibleWageClosure, model::Model{MobileLaborCES}, L_sum::Real, w::Real,
    cpi::Real) = L_sum - model.options.labor_bar

# ── Unit cost with the Cobb-Douglas limit guard (Stage 1.5, Milestone C) ──
# CES unit cost:  (A^(ϵ-1) · (fs·w^(1-ϵ) + (1-fs)·ip^(1-ϵ)))^(1/(1-ϵ)).
# At ϵ = 1 exactly the formula degenerates (x^Inf with x = 1 ± float error);
# the analytic Cobb-Douglas limit is  w^fs · ip^(1-fs) / A  and is used whenever
# |1 - ϵ| is below machine-relevant tolerance, keeping the CD special case
# sign-safe and continuous.
const _CD_LIMIT_TOL = 1e-9

function _ces_unit_cost(A_eff, fs, w, ip, ϵ)
    if abs(1 - ϵ) < _CD_LIMIT_TOL
        return (w .^ fs) .* (ip .^ (1 .- fs)) ./ A_eff
    end
    ((A_eff .^ (ϵ - 1)) .* (fs .* w .^ (1 - ϵ) .+ (1 .- fs) .* ip .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
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

const _ETA_SCALE_INDETERMINACY_TOL = 1.1e-6

# Project decision (ADR-0010): only the BF endpoints η = 0 (immobile baseline
# allocation) and η = 1 (fully cost-minimizing/mobile) are kept. The
# interpolated cases 0 < η < 1 were carried by the ad hoc allocation wedge
# (retired with this change); no intermediate allocation is supported.
function _checked_eta(η)
    isfinite(η) || throw(ArgumentError("η must be finite"))
    (η == 0.0 || η == 1.0) || throw(DomainError(η,
        "η must be 0 (immobile) or 1 (fully mobile); intermediate reallocation was retired (ADR-0010)"))
    Float64(η)
end

_positive_floor(x) = max.(x, eps(Float64))
_positive_floor(x::Real) = max(x, eps(Float64))

"""
	_intermediate_price(Ω_raw, p, θ)

Intermediate-goods price index. CES for θ ≠ 1; the exact Cobb-Douglas
limit ip_u = prod_j p_j^{Ω[u,j]} at θ = 1 (the bounded production-core
choice of Notebook 03b: with genuine sector self-loops -- up to Omega_ii
= 0.57 in the data -- the CES index at theta < 1 is self-referencing and
the zero-profit price system has exploding spiral branches; the
Cobb-Douglas index is log-linear with a unique positive root).
"""
function _intermediate_price(Ω_raw::AbstractMatrix, p::AbstractVector, θ::Real)
	isapprox(θ, 1.0; rtol = 0, atol = 1e-6) && return exp.(Ω_raw * log.(p))
	return (Ω_raw * p .^ (1 - θ)) .^ (1 / (1 - θ))
end

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

"""Sectoral labor allocation at the kept BF endpoints η ∈ {0, 1}.

η = 0 returns the baseline (immobile) allocation `data.labor_share`; η = 1
returns the cost-minimizing demand. The former interpolation (and the ad hoc
allocative-efficiency wedge that carried it) was retired with ADR-0010;
`_checked_eta` rejects intermediate values loudly.
"""
function _interpolated_labor(p, y, w, model::Model{MobileLaborCES})
    η = _checked_eta(model.options.elasticities.η)
    η == 1.0 && return _positive_floor(_cost_minimizing_labor(p, y, w, model))
    return _positive_floor(model.data.labor_share)
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
    _mobile_market_demand(model, p, y, w) -> NamedTuple

Goods-market blocks shared by `problem` (which enforces N-1 clearing equations
plus CPI = 1) and `market_clearing_residuals` (the full N-vector canary):

  intermediary_demand  sectoral intermediate demand
  total_final_demand   domestic final demand (household + additive + manna)
  c_dom                domestic household demand block (for the canary)
  additive             financed programme demand (for the canary)
  L_i                  sectoral labor allocation at the kept BF endpoints
  cost                 CES unit cost (no allocation wedge, ADR-0010)
  E                    household expenditure entering the CES demand

Budget-consistent CES demand over the household's (possibly financed)
expenditure. `fin = model.financing` (Foundation II):
  F1  composes E within its budget via preference weights (no additive
      demand; the normalizer keeps Σ p_i c_i = E exactly);
  F2  E = w·ΣL − T(p) with T = Σ p_i g_i the lump-sum tax;
  F3  E = w·ΣL, externally financed.
Compatibility (ADR-0005): the legacy unfinanced autonomous/investment manna
(A/G below) is RETAINED alongside the financed programme demand.
"""
function _mobile_market_demand(model::Model{MobileLaborCES}, p::AbstractVector,
        y::AbstractVector, w::Real)
    (; data, options, shocks) = model
    N = length(data.factor_share)
    (; consumption_share, Ω_raw, factor_share) = data
    (; θ, ϵ, σ) = options.elasticities
    intermediate_price = _intermediate_price(Ω_raw, p, θ)
    L_i = sectoral_labor_demand(p, y, w, model)
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    # NOTE: no positivity guard here — the residual function must tolerate the
    # solver's exploration of negative-income trial points (the legacy code
    # did). The E > 0 check belongs to the post-solve validation in the
    # notebooks (headline assertions verify E = w*L - T > 0 at equilibrium).
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(consumption_share .* ds_eff .* p .^ (1 - σ))
    # v3 open economy: the household consumes (1-s)E gross (saving s·E leaks);
    # only the DOMESTIC content circulates -- the import content of household,
    # government, investment and programme demand is supplied by the external
    # account. Exports are exogenous domestic sales (no margin).
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    cons_base = sum(data.labor_share)
    A = shocks.autonomous_demand .* data.consumption_share .* cons_base
    G = shocks.investment_shock .* data.consumption_share .* cons_base
    total_final_demand = c_dom .+
                         (1 .- data.import_margin) .* additive_demand(fin, N) .+
                         (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
                         data.exports_demand .+ A .+ G
    # Intermediate demand: the DOMESTIC bill A_u scaled by equilibrium output
    # (a_u = A_bill/λ_u; LaForge exact A-bill fix). The imported+taxed content
    # of the bill is an external-account leak (data.M_int), not domestic demand.
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (data.A_bill ./ data.λ) .* y))
    # With only the BF endpoints η ∈ {0, 1} (ADR-0010) the ad hoc B&F (2019)
    # "labor-reallocation wedge" is retired: it only existed to carry the
    # interpolated 0 < η < 1 cases and was never derived as a CES
    # allocative-loss coefficient.
    cost = _ces_unit_cost(shocks.supply_shock, factor_share, w, intermediate_price, ϵ)
    return (; intermediary_demand, total_final_demand, c_dom, additive = additive_demand(fin, N),
        L_i, cost, E)
end

"""
    problem(out, X, model::Model{MobileLaborCES})

The equilibrium system for the mobile-labor CES model.

Unknowns: X = [p(1:N); y(1:N); w]  — 2N+1 elements
Equations (2N+1):
  1. Zero-profit equations  (N):     p_i = cost_i(p, w)   for all i=1..N
  2. Market clearing        (N-1):   y_i = intermediary_demand_i + final_demand_i
                                      for i=1..N-1
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄
  4. Numeraire              (1):     CPI = 1  (Σ β_i · p_i^(1-σ))^(1/(1-σ) = 1)

Note (review finding 2.1, ADR-0010): the N-th market is NOT Walras-redundant
in the v3 open economy once imports leak — the p-weighted sum of the clearing
residuals equals the external-account imbalance. The system omits the N-th
clearing and defines the external balance residually: the omitted residual is
exposed by `market_clearing_residuals`, and the identity with
S − (I+X−M) is asserted by the acceptance tests/experiments. Enforcing all N
clearings instead (the all-N `644ba37` form, reverted on `revisefinal` by
`0f33ad6`) over-determines the open economy: the additive-demand cases then
have no exact root (measured floor 4.4e-4 on
the real 70-sector calibration, 9e-3 on the v3 contract fixture), because the
extra scale-invariant equation is not compatible with the accounting.

The CPI = 1 numeraire pins the price level (the mobile system is homogeneous
of degree 1 in (p, w) without it), mirroring the fixed-wage system
(`problem_fixed`), where w = 1 pins the scale.

Economic note: η selects the sectoral labor allocation at the two kept BF
endpoints — η = 0 keeps the baseline (immobile) allocation, η = 1 uses the
cost-minimizing (fully mobile) allocation (ADR-0010). It is not a
labor-supply elasticity.
"""
function problem(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options) = model
    N = length(data.factor_share)

    # FULL FORMULATION (2N+1 unknowns: p1..pN, y1..yN, w): ALL N zero-profit
    # conditions, N-1 clearing equations (the N-th market is the residual
    # external account; see the docstring and ADR-0010), the labour market and
    # the CPI = 1 numeraire.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)  # scalar wage, keep positive

    blocks = _mobile_market_demand(model, p, y, w)
    cpi = sum(data.consumption_share .* p .^ (1 - options.elasticities.σ))^(1 / (1 - options.elasticities.σ))

    # ── Equation 1: Zero-profit for ALL N sectors ──
    out[1:N] .= p .- blocks.cost

    # ── Equation 2: Market clearing for sectors 1..N-1 ──
    # The N-th clearing is NOT Walras-redundant once imports leak; it is the
    # residual external account and is checked by `market_clearing_residuals`
    # instead of being imposed here (ADR-0010).
    out[N+1:2N-1] .= y[1:N-1] .- blocks.intermediary_demand[1:N-1] .- blocks.total_final_demand[1:N-1]

    # ── Equation 3: Labour market (flexible-wage system: ALPHA / BETA) ──
    out[2N] = labor_market_residual(labor_closure(options), model, sum(blocks.L_i), w, cpi)

    # ── Equation 4: Numeraire constraint -- CPI = 1 ──
    out[2N+1] = cpi - 1.0

    nothing
end

"""
    market_clearing_residuals(model, X) -> Vector

Full N-vector of goods-market clearing residuals `y_i − inter_i − final_i` at
the mobile vector `X = [p; y; w]` (clamped exactly as in `problem`). The last
entry is the market that `problem` does NOT impose: it is the residual
external account (review finding 2.1, ADR-0010). At a solved system the other
N-1 entries are ~0, so `p ⋅ market_clearing_residuals` reduces to the omitted
market's residual and matches the external-account imbalance
`S − (I+X−M)` (asserted by the acceptance tests; the experiments record it as
the canary). Diagnostic only — never part of the solve.
"""
function market_clearing_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    length(X) == 2N + 1 || throw(DimensionMismatch(
        "mobile-labor canary expects a $(2N+1)-element vector [p; y; w]"))
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)
    blocks = _mobile_market_demand(model, p, y, w)
    return y .- blocks.intermediary_demand .- blocks.total_final_demand
end

"""
    external_balance_canary(model, X) -> NamedTuple

`S − (I + X − M) + T` at the mobile vector `X = [p; y; w]`, in value terms and
consistent with the model's demand blocks:
  S    = s·E,
  I+X  = p·(exo_demand + exports_demand),
  M    = import content of final demand: `m/(1-m)` on the domestic household
         block, `m` on the government/investment/programme injections, minus
         the full programme value under F3 (ExternalDebt finances it
         externally, so the inflow offsets the trade balance), plus the
         intermediate-import leak `M_int` (row 74, ADR-0012).
  T    = product taxes on intermediate use `T_int` (row 75, ADR-0013) — the
         third component of the purchaser-price intermediate bill, which the
         A-bill charges to no one: `(1−fs)·λ ≡ A_bill + M_int + T_int`.
At a mobile (η = 1) equilibrium the omitted N-th market residual
(`market_clearing_residuals`) equals this quantity; the acceptance tests
assert that identity (review finding 2.1, ADR-0010; the `T` term added by
ADR-0013 closes it to machine precision — measured 1.7e-16 on full-71).
At η = 0 the omitted market additionally reflects the fixed-allocation/
factor-market gap, and legacy manna (ADR-0005) has no modelled import content;
the identity is asserted only for the mobile, zero-manna case.
"""
function external_balance_canary(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    length(X) == 2N + 1 || throw(DimensionMismatch(
        "mobile-labor canary expects a $(2N+1)-element vector [p; y; w]"))
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)
    blocks = _mobile_market_demand(model, p, y, w)
    (; data) = model
    m = data.import_margin
    M_cons = dot(p .* (m ./ max.(1 .- m, eps(Float64))), blocks.c_dom)
    M_inj = dot(p .* m, blocks.additive .+ data.gov_demand .+ data.exo_demand)
    M_prog = model.financing isa ExternalDebt ? -dot(p, blocks.additive) : 0.0
    # Intermediate-bill leaks (A-bill fix): the two non-domestic components of
    # the purchaser-price intermediate bill, both scaling with sectoral output.
    # Row 74 (imported intermediates, ADR-0012) and row 75 (product taxes on
    # intermediate use, ADR-0013). Omitting row 75 leaves the identity short by
    # exactly that term (measured: -2.6e-2 on full-71 before the fix).
    M_intl = dot(p .* (data.M_int ./ data.λ), y)
    T_intl = dot(p .* (data.T_int ./ data.λ), y)
    S = data.saving_rate * blocks.E
    IX = dot(p, data.exo_demand .+ data.exports_demand)
    return (; S = S, IX = IX,
        M = M_cons + M_inj + M_prog + M_intl, T = T_intl,
        diff = S - (IX - (M_cons + M_inj + M_prog + M_intl + T_intl)))
end

"""Return the exact residual vector for either mobile-labor closure.

For the `:fixed` closure the canonical vector is FULL form (2N: p1..pN,
y1..yN; w = 1 pinned); a mobile (p, y, w) vector is also accepted and reduced
internally (the wage component is dropped).
"""
function equilibrium_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    fixed = labor_closure(model.options) isa FixedWageClosure
    Xv = collect(X)
    if fixed
        # Fixed: canonical FULL vector (2N: p1..pN, y1..yN; w = 1 pinned).
        if length(Xv) == 2N
            xr = Xv
        elseif length(Xv) == 2N + 1
            xr = Xv[1:2N]   # drop the wage component of a mobile vector
        else
            throw(DimensionMismatch("fixed closure expects a $(2N)-element vector"))
        end
        out = similar(xr)
        problem_fixed(out, xr, model)
    else
        # Mobile: canonical FULL vector (2N+1: p1..pN, y1..yN, w).
        length(Xv) == 2N + 1 || throw(DimensionMismatch(
            "mobile closure expects a $(2N+1)-element vector"))
        out = similar(Xv)
        problem(out, Xv, model)
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
w = 1.0 is hard-coded as the numeraire; unknowns are p(1:N) and y(1:N) only (2N).
Equations (2N):
  1. N zero-profit conditions:  p_i = cost_i(p, w=1.0) for all i
  2. N market-clearing:         y_i = intermed_i + final_i for all i
Total: 2N equations, 2N unknowns. There is no CPI pin: w = 1 is the numeraire.
Employment `L_i` is computed post-solve and is not constrained to `labor_bar`.
The sticky-wage counterfactual holds the production cost at its direct CES
form (the retired mobile interpolation wedge, ADR-0010, is not used here).
"""
function problem_fixed(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    N = length(data.factor_share)
    w = 1.0  # sticky wage

    # FIXED-WAGE FORMULATION (2N unknowns: p1..pN, y1..yN; w = 1 pinned as the
    # sticky-wage numeraire). ALL N zero-profit and ALL N clearing equations
    # are enforced -- with the v3 homogeneous budget there is no Walras
    # redundancy at fixed w, and no sector's zero-profit may be dropped.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])

    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω_raw, factor_share, labor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    # Intermediate goods price index
    intermediate_price = _intermediate_price(Ω_raw, p, θ)

    # CPI
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Sectoral labor demand at w=1.0
    L_i = sectoral_labor_demand(p, y, w, model)

    # Final demand (budget-consistent, financed — same hook as `problem`).
    fin = model.financing
    ds_eff = preference_weights(fin, demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    # NOTE: no positivity guard here — the residual function must tolerate the
    # solver's exploration of negative-income trial points (the legacy code
    # did). The E > 0 check belongs to the post-solve validation in the
    # notebooks (headline assertions verify E = w*L - T > 0 at equilibrium).
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(consumption_share .* ds_eff .* p .^ (1 - σ))
    # v3 open economy: the household consumes (1-s)E gross (saving s·E leaks);
    # only the DOMESTIC content circulates -- the import content of household,
    # government, investment and programme demand is supplied by the external
    # account. Exports are exogenous domestic sales (no margin). The leakages
    # s·E + M now balance the injections I + X -- the identity S = I + X - M
    # holds at every equilibrium (validated in the notebooks).
    c_dom = (1 .- data.saving_rate) .* (1 .- data.import_margin) .*
            (consumption_share .* ds_eff) .* E .* p .^ (-σ) ./ agg
    # Legacy compatibility (ADR-0005): the unfinanced autonomous/investment
    # manna is RETAINED alongside the financed programme demand. cbase2 retired
    # manna; the root kernel keeps it so the characterization goldens and
    # rerun_results.jl stay reproducible. With NoFinancing and zero v3 fields
    # the demand below is numerically identical to the pre-Phase-2 code.
    cons_base = sum(data.labor_share)
    A = shocks.autonomous_demand .* data.consumption_share .* cons_base
    G = shocks.investment_shock .* data.consumption_share .* cons_base
    total_final_demand = c_dom .+
                         (1 .- data.import_margin) .* additive_demand(fin, N) .+
                         (1 .- data.import_margin) .* (data.gov_demand .+ data.exo_demand) .+
                         data.exports_demand .+ A .+ G

    # Intermediary demand: DOMESTIC bill coefficient a_u = A_bill/λ_u (A-bill
    # fix); the imported+taxed content is an external-account leak.
    intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (data.A_bill ./ data.λ) .* y))

    # Direct CES cost at the sticky wage (CD-limit-safe).
    cost = _ces_unit_cost(supply_shock, factor_share, w, intermediate_price, ϵ)

    # 1. Zero-profit for ALL N sectors (p1 = 1 is NOT pinned here: the wage
    #    w = 1 is the numeraire of the sticky-wage regime, and every sector's
    #    zero-profit must hold)
    out[1:N] .= p .- cost

    # 2. Market clearing for ALL N sectors. At pinned w there is no Walras
    #    redundancy under the v3 homogeneous budget: the saving leak s·E
    #    balances the exogenous injections I + X (identity S = I + X - M).
    out[N+1:2N] .= y .- intermediary_demand .- total_final_demand

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
    # Scale determinacy of the fixed-wage η ≈ 1 system is a VERIFIED property of
    # the round-gain matrix, not a financing-type heuristic (ADR-0014). The
    # clearing block is y = G·y + const, with G's column sums equal to the
    # round-gain
    #     colsum_u = A_bill_u/λ_u + (1 − m_u)·(1 − s)·fs_u,
    # so (I − G) is singular — a unit root, hence a continuum of solutions —
    # exactly when the largest column sum reaches 1. On a CLOSED fixture
    # (m = s = 0, A_bill = (1−fs)λ) the column sums are exactly 1 and the guard
    # fires; on the OPEN A-bill calibration they are strictly below 1 (the
    # finiteness gate of `recalibrate_open`) and the system is determinate —
    # measured on full-71: σ_min/σ_max = 0.1586 at the F1 point, and three inits
    # (λ, 2λ, λ/2) converge to the same root (L = 0.9990271532561, agreeing to
    # 1e-15). The error keeps the legacy surfaces ("scale-indeterminate" and
    # "autonomous or investment") asserted by the contract tests.
    if isapprox(η, 1.0; rtol=0, atol=_ETA_SCALE_INDETERMINACY_TOL)
        colsums = data.A_bill ./ data.λ .+
            (1.0 .- data.import_margin) .* (1.0 - data.saving_rate) .* data.factor_share
        maximum(colsums) >= 1.0 - 1e-12 && throw(ArgumentError(
            "fixed-wage η=1 has a homogeneous, scale-indeterminate equilibrium " *
            "(max round-gain column sum = $(maximum(colsums)) ≥ 1, a unit root): " *
            "add autonomous or investment demand as an additive-demand anchor " *
            "(or a TaxFinanced / ExternalDebt programme bundle), or use another η"))
    end

    if init === nothing
        init = [ones(N); data.λ]                # p = 1, y = λ
    elseif length(init) == 2N + 1
        # Full (p, y, w) init from a mobile solution: drop the wage component
        init = init[1:2N]
    end

    # Avoid asking the nonlinear solver to differentiate an already exact
    # baseline, which some solver versions report as stalled.
    # Solver-internal failures (e.g. on non-finite trial points) are reported
    # as non-convergence to preserve the legacy "did not converge" surface.
    # (Length errors from the residual check propagate unchanged.)
    is_exact = maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
    x = if is_exact
        Float64.(init)
    else
        try
            ProbN = NonlinearSolve.NonlinearProblem(problem_fixed, init, model)
            res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=20000)
            # Quality gate = the ACTUAL residual, never the retcode. Bounded LM polish.
            x = res.u
            rmax = maximum(abs, equilibrium_residuals(model, x))
            for _ in 1:3
                rmax <= 1e-6 && break
                res = NonlinearSolve.solve(
                    NonlinearSolve.NonlinearProblem(problem_fixed, x, model),
                    NonlinearSolve.LevenbergMarquardt(); reltol=1e-8, abstol=1e-8, maxiters=20000)
                x = res.u
                rmax = maximum(abs, equilibrium_residuals(model, x))
            end
            if rmax > 1e-6
                error("MobileLaborCES._solve_fixed did not converge: retcode = $(res.retcode), max|resid| = $rmax")
            end
            x
        catch e
            e isa ErrorException && occursin("did not converge", sprint(showerror, e)) && rethrow()
            e isa DimensionMismatch && rethrow()
            error("MobileLaborCES._solve_fixed did not converge: $(sprint(showerror, e))")
        end
    end

    p = x[1:N]
    q = max.(x[N+1:2N], 0.0)
    w = 1.0  # sticky wage

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all 1.0 — sticky)
    wages = fill(w, N)

    # Consumption (budget-consistent, financed — must match `problem_fixed`)
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
    consumption = (1 .- data.saving_rate) .* (data.consumption_share .* ds_eff .* E .* p .^ (-σ)) ./ agg  # gross household consumption (saving sE leaks)

    # Real GDP: consumption Tornqvist (B&F metric). v3 base = the calibrated
    # baseline household block c0_gross (data.household_baseline) -- the exact
    # v3 baseline demand, so the index is 1 at the baseline by construction.
    # Intermediate homotopy/beta rungs may have E < 0 (negative consumption);
    # report NaN there -- those points are rejected by the callers anyway.
    base_consumption = data.household_baseline
    real_gdp_index = all(>=(0), consumption) ?
        tornqvist_quantity_index(p, consumption, ones(N), base_consumption) : NaN
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

    # ── Standard :mobile path (FULL 2N+1 system: p1..pN, y1..yN, w) ──
    if init === nothing
        # Default initialization: p=1, y=λ, wage=1
        init = [ones(N); data.λ; 1.0]
    end

    # As in the fixed-wage path, avoid asking the nonlinear solver to
    # differentiate an already exact baseline (some versions report this as
    # stalled). Solver-internal failures are reported as non-convergence to
    # preserve the legacy "did not converge" surface.
    # (Length errors from the residual check propagate unchanged.)
    is_exact = maximum(abs, equilibrium_residuals(model, init)) <= 1e-12
    x = if is_exact
        Float64.(init)
    else
        try
            ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
            res = NonlinearSolve.solve(ProbN, reltol=1e-6, abstol=1e-6, maxiters=20000)
            # Quality gate = the ACTUAL residual, never the retcode (NonlinearSolve
            # reports Stalled on slow final convergence). Bounded LM polish (up to
            # 3 attempts) from the last point; verify before accepting. Mobile gate:
            # 1e-5 -- at the near-singular labour-equation direction the FD-Newton
            # floor is ~3e-6 (sum L off by 0.0003 percent, economically nil); the
            # budget identities remain EXACT and the DELTA equivalence gate stays
            # at machine precision.
            x = res.u
            rmax = maximum(abs, equilibrium_residuals(model, x))
            for _ in 1:3
                rmax <= 1e-5 && break
                res = NonlinearSolve.solve(
                    NonlinearSolve.NonlinearProblem(problem, x, model),
                    NonlinearSolve.LevenbergMarquardt(); reltol=1e-8, abstol=1e-8, maxiters=20000)
                x = res.u
                rmax = maximum(abs, equilibrium_residuals(model, x))
            end
            if rmax > 1e-5
                error("MobileLaborCES.solve did not converge: retcode = $(res.retcode), max|resid| = $rmax")
            end
            x
        catch e
            e isa ErrorException && occursin("did not converge", sprint(showerror, e)) && rethrow()
            e isa DimensionMismatch && rethrow()
            error("MobileLaborCES.solve did not converge: $(sprint(showerror, e))")
        end
    end

    # Full unknowns -> Solution fields (clamp solver dust below zero; the
    # Tornqvist index requires nonnegative quantities)
    p = x[1:N]
    q = max.(x[N+1:2N], 0.0)
    w = x[2N+1]

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all equal to w — mobile labor)
    wages = fill(w, N)

    # Consumption — must match the budget-consistent, financed demand used
    # inside `problem` so the reported allocation is the actual equilibrium
    # allocation.
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = w * L_sum
    E = household_expenditure(fin, model, total_income, p, L_sum)
    agg = sum(data.consumption_share .* ds_eff .* p .^ (1 - σ))
    consumption = (1 .- data.saving_rate) .* (data.consumption_share .* ds_eff .* E .* p .^ (-σ)) ./ agg  # gross household consumption (saving sE leaks)

    # Real GDP: Tornqvist (Divisia) quantity index of FINAL CONSUMPTION (value-added)
    # — not gross output. In B&F (2019), real GDP is a Divisia index of real final
    # demand / value added. Using gross output confounds intermediate flows with
    # welfare-relevant final output and dilutes the reallocation bridge. The
    # consumption (final-demand) Tornqvist isolates the welfare-relevant change.
    # v3 base = the calibrated baseline household block (data.household_baseline),
    # so the index is 1 at the v3 baseline by construction.
    base_consumption = data.household_baseline
    real_gdp_index = all(>=(0), consumption) ?
        tornqvist_quantity_index(p, consumption, ones(N), base_consumption) : NaN
    nominal_gdp = w * sum(L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Utility: construct a model from data + elasticities
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=nothing, closure=:mobile, financing=nothing, eta_s=nothing)

Convenience constructor for a MobileLaborCES model.
`closure=:mobile` uses a flexible, market-clearing wage; `closure=:fixed` uses a
sticky wage of one and unconstrained employment demand; `closure=:beta` (or
`eta_s=<value>`, which forces `:beta`) uses the elastic-supply BETA closure.
`financing` selects an F1/F2/F3 financing closure (default `NoFinancing()`).
All existing positional/keyword forms keep working.
"""
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Real, Nothing}=nothing, closure=:mobile, financing::Union{AbstractFinancing, Nothing}=nothing, eta_s::Union{Real, Nothing}=nothing)
    # BETA: an explicit supply elasticity forces the :beta closure.
    if eta_s !== nothing
        closure in (:mobile, :beta) || throw(ArgumentError(
            "eta_s applies to the :beta closure (got closure = $closure)"))
        closure = :beta
    end
    el = eta_s === nothing ? MobileLaborCESElasticities(θ, ϵ, σ, η) :
                             MobileLaborCESElasticities(θ, ϵ, σ, η, Float64(eta_s))
    closure_symbol = _closure_symbol(closure)
    labor_bar !== nothing && closure_symbol == :fixed && throw(ArgumentError(
        "fixed closure treats employment as an outcome; labor_bar is not used and must not be supplied"))
    lb = labor_bar === nothing ? sum(data.labor_share) : Float64(labor_bar)
    fin = financing === nothing ? NoFinancing() : financing
    model = Model(data, shocks, MobileLaborCES(el, lb, closure_symbol), fin)
    return model
end

mobile_labor_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real;
    labor_bar=nothing, closure=:mobile, financing=nothing, eta_s=nothing) = mobile_labor_model(data, shocks, Float64(θ), Float64(ϵ), Float64(σ), Float64(η);
    labor_bar=labor_bar === nothing ? nothing : Float64(labor_bar),
    closure=_closure_symbol(closure),
    financing=financing,
    eta_s=eta_s)
