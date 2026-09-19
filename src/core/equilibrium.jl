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
	external_transfer::Float64
	model::Model
end

function Solution(prices_raw, quantities, wages, consumption, numeraire, real_gdp, nominal_gdp, model; external_transfer::Real = 0.0)
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
		Float64(external_transfer),
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
# Explicit external-account closure (ADR-0019, superseding ADR-0010's N-1
# decision): every regime enforces ALL N goods-market clearing equations. The
# mobile system gains one scalar unknown F (net external transfer, entering
# household expenditure after tax) alongside the CPI = 1 numeraire; the
# fixed-wage and η = 0 systems clear all N markets with F ≡ 0. The old all-N
# `644ba37` form (reverted on `revisefinal` by `0f33ad6`) had no exact root
# only because it lacked the external variable; with F the all-N system is
# well posed (baseline residual 4.4e-16). See `problem`. The allocation wedge
# stays retired with η ∈ {0,1}.

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
    # ADR-0022: the sectoral generalisation of BETA. `nothing` keeps every
    # single-wage path bit-identical; a length-N vector selects the N sectoral
    # labour markets (the 3N+1 system `problem_sectoral`, whose eta_s,i = 0
    # corner is exactly the ADR-0020 option C endpoint).
    eta_s_vec::Union{Nothing,Vector{Float64}}
end
# Parent-compatible 4-arg constructor (no elastic labour supply).
MobileLaborCESElasticities(θ::Real, ϵ::Real, σ::Real, η::Real) =
    MobileLaborCESElasticities(Float64(θ), Float64(ϵ), Float64(σ), Float64(η), 0.0, nothing)
# Scalar BETA (5-arg): the single-wage system, unchanged.
MobileLaborCESElasticities(θ::Real, ϵ::Real, σ::Real, η::Real, eta_s::Real) =
    MobileLaborCESElasticities(Float64(θ), Float64(ϵ), Float64(σ), Float64(η),
        Float64(eta_s), nothing)
# Sectoral BETA (ADR-0022): the length-N elasticity vector.
MobileLaborCESElasticities(θ::Real, ϵ::Real, σ::Real, η::Real, eta_s::Real,
        eta_s_vec::AbstractVector) =
    MobileLaborCESElasticities(Float64(θ), Float64(ϵ), Float64(σ), Float64(η),
        Float64(eta_s), Float64.(collect(eta_s_vec)))

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
    options.closure == :beta   ? (options.elasticities.eta_s_vec === nothing ?
        ElasticLaborClosure(options.elasticities.eta_s) :
        SectoralElasticLaborClosure(options.elasticities.eta_s_vec)) :
                                 FixedWageClosure()

# ── Labor-market equation hook ──
# The total-labor-market residual of the flexible-wage systems (problem(), 2N+2).
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

# ADR-0017: the fixed-wage η ≈ 1 admissibility guard assesses the ACTUAL
# clearing matrix G of `problem_fixed` (`y = G·y + c` at the baseline
# reference prices p = ones). The strict column-sum bound is only a
# SUFFICIENT contraction shortcut; when it fails, determinacy is decided by
# the rank of (I − G) and the sign of the unique candidate solution.
const _FIXED_CLEARING_COLSUM_TOL = 1e-12
const _FIXED_CLEARING_RANK_TOL = 1e-10
const _FIXED_CLEARING_POSITIVITY_TOL = 1e-12

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

"""Cost-minimizing labor demand, evaluated in log space for stable η extrapolation.

`w` is a scalar (one economy-wide wage, every η = 1 regime) or a vector (the
sectoral wages of the η = 0 endpoint, ADR-0020 option C). The expression is
elementwise in `w`, so the scalar case is unchanged to the last bit.
"""
function _cost_minimizing_labor(p, y, w, model::Model{MobileLaborCES})
    (; data, options, shocks) = model
    (; ϵ) = options.elasticities
    (; factor_share) = data
    p, y, w, A, α = _positive_floor.((p, y, w, shocks.supply_shock, factor_share))
    log_demand = ϵ .* (log.(p) .+ ((ϵ - 1) / ϵ) .* log.(A) .+ (1 / ϵ) .* log.(α) .+
        (1 / ϵ) .* log.(y) .- log.(w))
    exp.(clamp.(log_demand, log(floatmin(Float64)), log(floatmax(Float64))))
end

"""
	_wage_bill(w, L)

Household wage income. The scalar-wage path keeps its exact arithmetic
(`w * sum(L)`); the sectoral-wage path (ADR-0020, η = 0) sums `w_i * L_i`.
"""
_wage_bill(w::Real, L) = w * sum(L)
_wage_bill(w::AbstractVector, L) = sum(w .* L)

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
# The equilibrium problem (2N+2 equations, 2N+2 unknowns)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _mobile_market_demand(model, p, y, w; external_transfer = 0.0) -> NamedTuple

Goods-market blocks shared by `problem` (all N clearing equations plus the
labour market and CPI = 1), `problem_fixed` (F ≡ 0), `market_clearing_residuals`
and the canary/diagnostics:

  intermediary_demand  sectoral intermediate demand
  total_final_demand   domestic final demand (household + additive + manna)
  c_dom                domestic household demand block (for the canary)
  additive             financed programme demand (for the canary)
  L_i                  sectoral labor allocation at the kept BF endpoints
  cost                 CES unit cost (no allocation wedge, ADR-0010)
  E                    household expenditure entering the CES demand
  (after-tax wage income plus the external transfer F)

Budget-consistent CES demand over the household's (possibly financed)
expenditure. `fin = model.financing` (Foundation II):
  F1  composes E within its budget via preference weights (no additive
      demand; the normalizer keeps Σ p_i c_i = E exactly);
  F2  E = w·ΣL − T(p) + F with T = Σ p_i g_i the lump-sum tax (F after tax);
  F3  E = w·ΣL + F, externally financed.
Compatibility (ADR-0005): the legacy unfinanced autonomous/investment manna
(A/G below) is RETAINED alongside the financed programme demand.

`w` is a scalar in every η = 1 regime (one economy-wide wage) and the sectoral
wage vector at the η = 0 endpoint (`problem_sectoral`, ADR-0020 option C):
there `L_i` is the frozen `data.labor_share` (returned by
`sectoral_labor_demand` at η = 0, which ignores `w`) and household wage income
is `sum_i w_i L_i` (`_wage_bill`). The scalar-wage arithmetic is unchanged.
"""
function _mobile_market_demand(model::Model{MobileLaborCES}, p::AbstractVector,
        y::AbstractVector, w::Union{Real,AbstractVector}; external_transfer::Real = 0.0)
    (; data, options, shocks) = model
    N = length(data.factor_share)
    (; consumption_share, Ω_raw, factor_share) = data
    (; θ, ϵ, σ) = options.elasticities
    intermediate_price = _intermediate_price(Ω_raw, p, θ)
    L_i = sectoral_labor_demand(p, y, w, model)
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = _wage_bill(w, L_i)
    # NOTE: no positivity guard here — the residual function must tolerate the
    # solver's exploration of negative-income trial points (the legacy code
    # did). The E > 0 check belongs to the post-solve validation in the
    # notebooks (headline assertions verify E = w*L - T > 0 at equilibrium).
    E = household_expenditure(fin, model, total_income, p, L_sum; external_transfer = external_transfer)
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

Unknowns: X = [p(1:N); y(1:N); w; F]  — 2N+2 elements, where F is the net
external transfer (net external borrowing) entering household expenditure
after tax: E = (1-τ)·w·ΣL + F (ADR-0019, decision D1).
Equations (2N+2):
  1. Zero-profit equations  (N):     p_i = cost_i(p, w)   for all i=1..N
  2. Market clearing        (N):     y_i = intermediary_demand_i + final_demand_i
                                      for all i=1..N
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄ (ALPHA), or the BETA
                                      elastic-supply curve; at η = 0 the PIN
                                      F = 0 (see below)
  4. Numeraire              (1):     CPI = 1  (Σ β_i · p_i^(1-σ))^(1/(1-σ) = 1)

At every all-N solution the exact external-account identity holds:
S + T_int + M − (I+X) = F + B_gov, where B_gov = Σ p_i g_i under ExternalDebt
(F3) financing and 0 otherwise, S = s·E with E including F, and M is the full
import content (final margin plus the intermediate-import leak M_int, row 74).
`external_balance_canary` asserts it (the acceptance test, at machine precision).
With nonzero legacy manna (the ADR-0005 compatibility path) the unfinanced
demand is not booked in B_gov: at an all-N η = 1 solution the gap reads exactly
`p·(A+G)` and F absorbs it (`F = −p·(A+G)` at zero programme). The matrix
designs pass zero manna, where the identity above holds as written.

Note (ADR-0019, superseding ADR-0010): no market is omitted. The old note —
that the all-N `644ba37` form, reverted on `revisefinal` by `0f33ad6`,
over-determines the open economy and has no exact root for additive-demand
cases — described a system that lacked the external variable; with F the
all-N formulation is well posed (measured baseline residual 4.4e-16 with a
well-conditioned Jacobian).

At η = 0 (the BF immobile endpoint) this system is NOT used: `sectoral_labor_demand`
returns the constant `data.labor_share`, so the labour residual carries no
information, and the endpoint is formulated as the sectoral-wage system
`problem_sectoral` (ADR-0020 option C), which replaces the retired `F = 0` pin
of ADR-0019 D2. Calling `problem` on an η = 0 model throws.

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

    # FULL FORMULATION (2N+2 unknowns: p1..pN, y1..yN, w, F): ALL N zero-profit
    # conditions, ALL N clearing equations, the labour market (or the F = 0 pin
    # at eta = 0, where the labour row is identically satisfied; decision D2)
    # and the CPI = 1 numeraire.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = max(X[2N+1], 1e-10)  # scalar wage, keep positive
    F = X[2N+2]  # net external transfer (signed: no floor)

    blocks = _mobile_market_demand(model, p, y, w; external_transfer = F)
    cpi = sum(data.consumption_share .* p .^ (1 - options.elasticities.σ))^(1 / (1 - options.elasticities.σ))

    # ── Equation 1: Zero-profit for ALL N sectors ──
    out[1:N] .= p .- blocks.cost

    # ── Equation 2: Market clearing for ALL N sectors ──
    out[N+1:2N] .= y .- blocks.intermediary_demand .- blocks.total_final_demand

    # ── Equation 3: Labour market (flexible-wage system: ALPHA / BETA) ──
    # The eta = 0 endpoint has its own formulation, `problem_sectoral`
    # (ADR-0020 option C); this system is the eta = 1 mobile one only.
    options.elasticities.η == 0.0 && throw(ArgumentError(
        "the eta = 0 endpoint is the sectoral-wage system (ADR-0020, option C); call problem_sectoral"))
    out[2N+1] = labor_market_residual(labor_closure(options), model, sum(blocks.L_i), w, cpi)

    # ── Equation 4: Numeraire constraint -- CPI = 1 ──
    out[2N+2] = cpi - 1.0

    nothing
end

"""
    problem_sectoral(out, X, model::Model{MobileLaborCES})

The η = 0 (BF) endpoint under ADR-0020 option C: **sectoral wages with the
frozen allocation**. Labour cannot reallocate, so the sectoral allocation stays
at `data.labor_share` and each sector's wage is set by that sector's own
marginal product at the frozen allocation.

Unknowns: `X = [p(1:N); y(1:N); w(1:N); F]` — 3N + 1.
Equations (3N + 1):
  1. zero-profit per sector (N):      `p_i = cost_i(p, w_i)`
  2. sectoral FOC at the frozen allocation (N):
                                      `log L^cm_i(p_i, y_i, w_i) = log labor_share_i`
  3. all-N clearing (N):              `y_i = inter_i + final_i`, household wage
                                      income `sum_i w_i L_i`,
                                      `E = (1 - tau) sum_i w_i L_i + F`
  4. numeraire (1):                   CPI = 1

`F` is the free scalar the equation count requires. The block of equations 1-3
is homogeneous of degree 1 in `(p, w, F)`, so its 3N equations determine 3N - 1
effective unknowns; replacing the mobile system's single aggregate labour
equation by N sectoral conditions adds N - 1 equations, so the block carries one
equation more than it has directions to pin and the demand side needs one free
scalar to be consistent with the supply side. `F` is that scalar (it enters `E`
after tax, exactly as in ADR-0019). Dropping one clearing equation instead is
the ADR-0010 shortcut that ADR-0019 retired, so it is not available.

Properties (both measured on the full-71 A-bill calibration, 2026-09-18):
at a solution the frozen allocation is cost-minimizing at these wages, so the
external-account identity gap `S + T_int + M - (I+X) - (F + B_gov)`, which is
exactly `-(sum_i w_i L^cm_i - sum_i w_i L_i)`, vanishes (the retired pin left up
to -0.79 percent of GDP open); and financing neutrality extends to η = 0
(`F_F3 = F_F2 - B_gov`, identical real allocation and booked position).
"""
function problem_sectoral(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data, options) = model
    N = length(data.factor_share)

    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])
    w = _positive_floor(X[2N+1:3N])
    F = X[3N+1]

    blocks = _mobile_market_demand(model, p, y, w; external_transfer = F)
    cpi = sum(data.consumption_share .* p .^ (1 - options.elasticities.σ))^(1 / (1 - options.elasticities.σ))

    # ── Equation 1: zero-profit for ALL N sectors, sector-specific wage ──
    out[1:N] .= p .- blocks.cost

    # ── Equation 2: sectoral labour-market conditions ──
    # The COST-MINIMIZING demand at this sector's wage must equal the sectoral
    # supply. `blocks.L_i` is NOT used here: at eta = 0
    # `sectoral_labor_demand` returns the frozen `labor_share` itself, which
    # would make this row identically zero and the system degenerate (the
    # external-account identity gate catches that, since the gap is exactly
    # the value of this row).
    #   ADR-0020 option C  : L^cm_i = Lbar_i                  (eta_s,i = 0)
    #   ADR-0022           : L^cm_i = Lbar_i * (w_i/CPI)^{eta_s,i}
    # The eta_s,i = 0 branch is kept arithmetically separate so that every
    # existing eta = 0 cell stays bit-identical.
    Lcm = _cost_minimizing_labor(p, y, w, model)
    esv = options.elasticities.eta_s_vec
    if esv === nothing
        out[N+1:2N] .= log.(Lcm) .- log.(_positive_floor(data.labor_share))
    else
        length(esv) == N || throw(DimensionMismatch(
            "the sectoral elasticity vector must have $N entries (got $(length(esv)))"))
        out[N+1:2N] .= log.(Lcm) .- log.(_positive_floor(data.labor_share)) .-
                       esv .* log.(w ./ cpi)
    end

    # ── Equation 3: market clearing for ALL N sectors ──
    out[2N+1:3N] .= y .- blocks.intermediary_demand .- blocks.total_final_demand

    # ── Equation 4: Numeraire constraint -- CPI = 1 ──
    out[3N+1] = cpi - 1.0

    nothing
end

"""
    market_clearing_residuals(model, X) -> Vector

Full N-vector of goods-market clearing residuals `y_i − inter_i − final_i`,
with the external transfer threaded into the demand hook. X is [p; y] (2N,
w = 1, F = 0), [p; y; w] (2N+1 legacy, F = 0) or [p; y; w; F] (2N+2),
clamped exactly as in `problem`. At a solved all-N system every entry is ~0
(ADR-0019: no market is omitted any more). Diagnostic only — never part of
the solve.
"""
function market_clearing_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    Xv = collect(X)
    if length(Xv) == 2N
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = 1.0
        F = 0.0
    elseif length(Xv) == 2N + 1
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = max(Xv[2N+1], 1e-10)
        F = 0.0
    elseif length(Xv) == 2N + 2
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = max(Xv[2N+1], 1e-10)
        F = Xv[2N+2]
    elseif length(Xv) == 3N + 1
        # eta = 0 sectoral-wage vector [p; y; w(1:N); F] (ADR-0020 option C)
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = _positive_floor(Xv[2N+1:3N])
        F = Xv[3N+1]
    else
        throw(DimensionMismatch(
            "mobile-labor canary expects a $(2N)-, $(2N+1)-, $(2N+2)- or $(3N+1)-element vector [p; y(; w(; F))]"))
    end
    blocks = _mobile_market_demand(model, p, y, w; external_transfer = F)
    return y .- blocks.intermediary_demand .- blocks.total_final_demand
end

"""
    external_balance_canary(model, X) -> NamedTuple

External-account identity gap at [p; y] (2N, w = 1, F = 0), [p; y; w]
(2N+1 legacy, F = 0) or [p; y; w; F] (2N+2), in value terms and consistent
with the model's demand blocks (the transfer F is threaded into the hook,
so S = s·E includes it):
  S    = s·E,
  I+X  = p·(exo_demand + exports_demand),
  M    = full import content: `m/(1-m)` on the domestic household block, `m`
         on the government/investment/programme injections, plus the
         intermediate-import leak `M_int` (row 74, ADR-0012). The programme's
         own import content legitimately belongs to M (it is no longer netted
         out); the external financing of the programme is booked in B_gov.
  T    = product taxes on intermediate use `T_int` (row 75, ADR-0013) — the
         third component of the purchaser-price intermediate bill, which the
         A-bill charges to no one: `(1−fs)·λ ≡ A_bill + M_int + T_int`.
         `M_int` and `T_int` are constant shares of the SAME per-user CES
         intermediate bundle the domestic bill `A_bill` is charged from, so
         both leaks are valued with the CES bill factor
         `k_u = p_u^ϵ · a_u^(ϵ−1) · P_u^(1-ϵ)` (ADR-0016),
  B_gov = Σ p_i g_i under ExternalDebt (F3) financing — the programme value —
         and 0 otherwise,
  diff = S + T + M − (I+X) − (F + B_gov): the external-account identity gap.
         With zero legacy manna (the matrix designs; ADR-0005 path) it is
         ≈ 0 at every η = 1 solution (mobile and fixed-wage, where the
         cost-minimizing allocation closes the account); at BF η = 0 it
         carries the fixed-allocation factor-market gap instead (measured
         4e-4..8e-3 on full-71). With nonzero legacy manna the unfinanced
         demand p·(A+G) is not booked in B_gov, so the gap reads exactly
         `p·(A+G)` at an all-N solution and F absorbs it (measured: F =
         −p·(A+G), diff = +p·(A+G) on a 3-sector fixture). `financing =
         F + B_gov` is the booked external position. The companion exact
         identity is `gdp_components(...).wedge = −diff`, on and off
         equilibrium.
The old reading (the omitted-market residual equals the imbalance) is
superseded by ADR-0019: every market clears, and the identity closes through
the booked external position instead.
"""
function external_balance_canary(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    Xv = collect(X)
    if length(Xv) == 2N
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = 1.0
        F = 0.0
    elseif length(Xv) == 2N + 1
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = max(Xv[2N+1], 1e-10)
        F = 0.0
    elseif length(Xv) == 2N + 2
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = max(Xv[2N+1], 1e-10)
        F = Xv[2N+2]
    elseif length(Xv) == 3N + 1
        # eta = 0 sectoral-wage vector [p; y; w(1:N); F] (ADR-0020 option C)
        p = _positive_floor(Xv[1:N])
        y = _positive_floor(Xv[N+1:2N])
        w = _positive_floor(Xv[2N+1:3N])
        F = Xv[3N+1]
    else
        throw(DimensionMismatch(
            "mobile-labor canary expects a $(2N)-, $(2N+1)-, $(2N+2)- or $(3N+1)-element vector [p; y(; w(; F))]"))
    end
    blocks = _mobile_market_demand(model, p, y, w; external_transfer = F)
    (; data, options, shocks) = model
    (; θ, ϵ) = options.elasticities
    m = data.import_margin
    M_cons = dot(p .* (m ./ max.(1 .- m, eps(Float64))), blocks.c_dom)
    M_inj = dot(p .* m, blocks.additive .+ data.gov_demand .+ data.exo_demand)
    # Intermediate-bill leaks (A-bill fix): the two non-domestic components of
    # the purchaser-price intermediate bill. They are constant shares of the
    # SAME per-user CES intermediate bundle the domestic bill A_bill is charged
    # from, so both carry the CES bill factor k_u = p_u^ϵ · a_u^(ϵ−1) · P_u^(1-ϵ)
    # (ADR-0016), where a is the supply shock and P = _intermediate_price.
    # Row 74 (imported intermediates, ADR-0012) and row 75 (product taxes on
    # intermediate use, ADR-0013). Omitting row 75 leaves the identity short by
    # exactly that term (measured: -2.6e-2 on full-71 before the fix).
    k_bill = p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* _intermediate_price(data.Ω_raw, p, θ) .^ (1 - ϵ)
    M_intl = dot(k_bill .* (data.M_int ./ data.λ), y)
    T_intl = dot(k_bill .* (data.T_int ./ data.λ), y)
    S = data.saving_rate * blocks.E
    IX = dot(p, data.exo_demand .+ data.exports_demand)
    B_gov = model.financing isa ExternalDebt ? dot(p, blocks.additive) : 0.0
    M = M_cons + M_inj + M_intl
    return (; S = S, IX = IX,
        M = M, T = T_intl,
        diff = S + T_intl + M - IX - (F + B_gov),
        external_transfer = F, programme_financing = B_gov, financing = F + B_gov)
end

"""Return the exact residual vector for either mobile-labor closure.

For the `:fixed` closure the canonical vector is FULL form (2N: p1..pN,
y1..yN; w = 1 pinned, F ≡ 0); a mobile (p, y, w) vector (2N+1) or a
(p, y, w, F) vector (2N+2) is also accepted and reduced internally (the wage
and transfer components are dropped).
For the mobile closures at η = 1 the canonical vector is (2N+2: p1..pN, y1..yN,
w, F); a legacy (p, y, w) vector (2N+1) is accepted with F = 0 for backward
compatibility with stored v1-v5 solutions/diagnostics.
At η = 0 the canonical vector is the sectoral-wage form (3N+1: p1..pN, y1..yN,
w1..wN, F) of `problem_sectoral` (ADR-0020 option C); a legacy 2N+2 vector
(scalar wage, the retired pin form) is expanded by replicating the wage.
"""
function equilibrium_residuals(model::Model{MobileLaborCES}, X::AbstractVector)
    N = length(model.data.factor_share)
    fixed = labor_closure(model.options) isa FixedWageClosure
    Xv = collect(X)
    if fixed
        # Fixed: canonical FULL vector (2N: p1..pN, y1..yN; w = 1 pinned, F ≡ 0).
        if length(Xv) == 2N
            xr = Xv
        elseif length(Xv) == 2N + 1
            xr = Xv[1:2N]   # drop the wage component of a mobile vector
        elseif length(Xv) == 2N + 2
            xr = Xv[1:2N]   # drop the wage and transfer components
        else
            throw(DimensionMismatch("fixed closure expects a $(2N)-element vector"))
        end
        out = similar(xr)
        problem_fixed(out, xr, model)
    elseif model.options.elasticities.η == 0.0 ||
           model.options.elasticities.eta_s_vec !== nothing
        # Sectoral system (ADR-0020 option C, ADR-0022): canonical 3N+1 vector
        # [p; y; w(1:N); F]. A legacy 2N+2 vector (scalar wage) is expanded by
        # replicating the wage.
        length(Xv) == 2N + 2 && (Xv = [Xv[1:2N]; fill(Xv[2N+1], N); Xv[2N+2]])
        length(Xv) == 3N + 1 || throw(DimensionMismatch(
            "the sectoral-wage closure expects a $(3N+1)-element vector [p; y; w; F]"))
        out = similar(Xv)
        problem_sectoral(out, Xv, model)
    else
        # Mobile (eta = 1): canonical FULL vector (2N+2: p1..pN, y1..yN, w, F).
        if length(Xv) == 2N + 1
            Xv = [Xv; 0.0]  # legacy vector: F = 0
        end
        length(Xv) == 2N + 2 || throw(DimensionMismatch(
            "mobile closure expects a $(2N+2)-element vector"))
        out = similar(Xv)
        problem(out, Xv, model)
    end
    out
end

function _equilibrium_residuals(model::Model{MobileLaborCES}, sol::Solution)
    X = labor_closure(model) isa FixedWageClosure ?
        [sol.prices_raw; sol.quantities] :
        (model.options.elasticities.η == 0.0 ||
         model.options.elasticities.eta_s_vec !== nothing) ?
        [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer] :
        [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
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
The fixed regime clears with F ≡ 0 (ADR-0019): there is no external-transfer
unknown and the demand hook below runs at its F = 0 default.
Employment `L_i` is computed post-solve and is not constrained to `labor_bar`.
The sticky-wage counterfactual holds the production cost at its direct CES
form (the retired mobile interpolation wedge, ADR-0010, is not used here).
The demand block is evaluated by `_mobile_market_demand` (w = 1), the
single demand kernel shared with the mobile system and the ADR-0017
admissibility guard, so the guard's clearing matrix cannot drift from
these equations.
"""
function problem_fixed(out::Vector, X::Vector, model::Model{MobileLaborCES})
    (; data) = model
    N = length(data.factor_share)

    # FIXED-WAGE FORMULATION (2N unknowns: p1..pN, y1..yN; w = 1 pinned as the
    # sticky-wage numeraire). ALL N zero-profit and ALL N clearing equations
    # are enforced -- with the v3 homogeneous budget there is no Walras
    # redundancy at fixed w, and no sector's zero-profit may be dropped.
    # The demand block is evaluated by `_mobile_market_demand` (w = 1), the
    # single demand kernel shared with the mobile system and the ADR-0017
    # admissibility guard, so the guard's clearing matrix cannot drift from
    # these equations.
    p = _positive_floor(X[1:N])
    y = _positive_floor(X[N+1:2N])

    blocks = _mobile_market_demand(model, p, y, 1.0)
    out[1:N] .= p .- blocks.cost
    out[N+1:2N] .= y .- blocks.intermediary_demand .- blocks.total_final_demand

    nothing
end

"""
    _fixed_clearing_affine(model, p) -> (G, c)

Affine representation `intermediary_demand + total_final_demand = G·y + c` of
the fixed-wage clearing block at prices `p` (w = 1). At fixed prices the block
is exactly affine in `y` — intermediate demand is linear in `y`, and household
expenditure/consumption are affine in `Σ L_i` — so `G` is extracted exactly
(to round-off) from unit differences at `y = λ` through the same demand hook
(`_mobile_market_demand`) that `problem_fixed` evaluates. The ADR-0017
admissibility guard calls this at the calibration baseline `p = ones(N)`;
callers must pass positive prices and positive reference output.
"""
function _fixed_clearing_affine(model::Model{MobileLaborCES}, p::AbstractVector)
    N = length(model.data.factor_share)
    y_ref = model.data.λ
    clearing(y) = begin
        blocks = _mobile_market_demand(model, p, y, 1.0)
        blocks.intermediary_demand .+ blocks.total_final_demand
    end
    b_ref = clearing(y_ref)
    G = Matrix{Float64}(undef, N, N)
    y_trial = copy(y_ref)
    for u in 1:N
        y_trial[u] = y_ref[u] + 1.0
        G[:, u] .= clearing(y_trial) .- b_ref
        y_trial[u] = y_ref[u]
    end
    return G, b_ref .- G * y_ref
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
    # Scale determinacy of the fixed-wage η ≈ 1 system is assessed on the
    # ACTUAL clearing matrix G of `problem_fixed` (`y = G·y + c`), built by
    # `_fixed_clearing_affine` from the same demand hook at the calibration
    # baseline prices (p = ones). ADR-0017 supersedes the ADR-0014 closed-form
    # round-gain trigger, which was wrong twice: a maximum column sum of 1 does
    # not imply a unit root, and the closed-form sums omitted the household
    # expenditure composition (import margins enter at the spending sector, so
    # the actual column sum is (1−s)·fs_u·Σ_i((1−m_i)·b_i), not
    # (1−m_u)(1−s)fs_u). The assessment distinguishes three properties:
    #   1. CONTRACTION: max(colsum(G)) < 1 is a SUFFICIENT shortcut for
    #      ρ(G) < 1 (G ≥ 0 at the reference) → (I − G) nonsingular, unique
    #      fixed point → admit without further checks;
    #   2. NONSINGULARITY: otherwise (I − G) is tested; a singular (I − G) is
    #      a unit root — a continuum of solutions — and the cell is rejected
    #      (the legacy "scale-indeterminate" surface);
    #   3. POSITIVE-SOLUTION EXISTENCE: a nonsingular (I − G) has the unique
    #      candidate y* = (I − G)⁻¹c at the reference point; it is admitted
    #      only if y* is positive at the model's resolution (an all-zero or
    #      sign-crossing candidate is admitted by nothing: the solver would
    #      return a floored, path-dependent point).
    if isapprox(η, 1.0; rtol=0, atol=_ETA_SCALE_INDETERMINACY_TOL)
        # Assess the η = 1 endpoint system: near-one η is snapped to 1.0 so
        # the shared demand hook — which admits only the kept BF endpoints
        # (ADR-0010) — evaluates the cost-minimizing allocation the guard
        # adjudicates. This preserves the legacy surface: near-one η on a
        # closed fixture throws the scale-indeterminacy ArgumentError before
        # any endpoint validation runs.
        el = options.elasticities
        model1 = η == 1.0 ? model : Model(model.data, model.shocks,
            MobileLaborCES(MobileLaborCESElasticities(el.θ, el.ϵ, el.σ, 1.0, el.eta_s),
                options.labor_bar, options.closure),
            model.financing)
        G, c = _fixed_clearing_affine(model1, ones(N))
        colsums = vec(sum(G; dims = 1))
        if maximum(colsums) >= 1.0 - _FIXED_CLEARING_COLSUM_TOL
            F = I - G
            sv = svdvals(F)
            if minimum(sv) <= _FIXED_CLEARING_RANK_TOL * maximum(sv)
                throw(ArgumentError(
                    "fixed-wage η=1 is scale-indeterminate: the actual clearing " *
                    "matrix (I − G) is singular at the baseline prices " *
                    "(σ_min/σ_max = $(minimum(sv) / maximum(sv)) ≤ " *
                    "$(_FIXED_CLEARING_RANK_TOL)), so (I − G) has a unit root " *
                    "and the equilibrium set contains a continuum. Additive " *
                    "demand — autonomous or investment manna, or a " *
                    "TaxFinanced / ExternalDebt programme bundle — is a constant " *
                    "and cannot remove the unit root; give the calibration a " *
                    "genuine leakage (positive saving or import margin), or use " *
                    "another η"))
            end
            y_star = F \ c
            scale = max(1.0, maximum(abs, y_star))
            if minimum(y_star) <= _FIXED_CLEARING_POSITIVITY_TOL * scale
                throw(ArgumentError(
                    "fixed-wage η=1 is determinate but has no positive " *
                    "equilibrium at the baseline prices: the unique clearing " *
                    "solution y = (I − G)⁻¹c is not positive " *
                    "(min y = $(minimum(y_star)) ≤ " *
                    "$(_FIXED_CLEARING_POSITIVITY_TOL)·|y|_max). Check the " *
                    "calibration / financing, or use another η"))
            end
        end
    end

    if init === nothing
        init = [ones(N); data.λ]                # p = 1, y = λ
    elseif length(init) == 2N + 1
        # Full (p, y, w) init from a mobile solution: drop the wage component
        init = init[1:2N]
    elseif length(init) == 2N + 2
        # Full (p, y, w, F) init from a mobile solution: drop wage and transfer
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
            # ADR-0015: polish to ~1e-10 even when the primary solve already sits
            # inside the 1e-6 acceptance gate. Cells that stop at the gate have
            # TOLERANCE-DEPENDENT metrics (GAMMA-F2 moved 8.3e-7 in real_gdp_rel
            # under a 1e-16 perturbation at resid 4.85e-7). The acceptance gate
            # below is unchanged, and the polish is monotone: a step that does
            # not improve the residual is discarded, so polishing can never turn
            # a passing cell into a failing one.
            for _ in 1:4
                rmax <= 1e-10 && break
                res = NonlinearSolve.solve(
                    NonlinearSolve.NonlinearProblem(problem_fixed, x, model),
                    NonlinearSolve.LevenbergMarquardt(); reltol=1e-12, abstol=1e-12, maxiters=20000)
                x_new = res.u
                r_new = maximum(abs, equilibrium_residuals(model, x_new))
                r_new < rmax || break
                x, rmax = x_new, r_new
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
    E = household_expenditure(fin, model, total_income, p, L_sum; external_transfer = 0.0)
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

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model; external_transfer = 0.0)
end


"""
    solve(model::Model{MobileLaborCES}; init)

Solve the mobile-labor CES model. Returns a `Solution` with the equilibrium
prices, quantities, wage (as a vector of the same wage in all sectors),
external transfer F, consumption, and GDP measures.

The initial guess defaults to baseline prices=1, quantities=λ, wage=1 and
F=0. A legacy 2N+1 init (without F) is accepted with F = 0.
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

    # ── Standard :mobile path ──
    # eta = 1: FULL 2N+2 system (p1..pN, y1..yN, w, F).
    # eta = 0: the sectoral-wage endpoint (ADR-0020 option C), canonical 3N+1
    # vector [p; y; w(1:N); F]. A warm start given in the 2N+2 mobile form
    # (scalar wage) is expanded by replicating the wage.
    η0 = options.elasticities.η == 0.0
    # ADR-0022: a sectoral elasticity vector selects the N-market system, which
    # is the eta = 0 formulation generalised by the real-wage supply term.
    esv = options.elasticities.eta_s_vec
    if esv !== nothing && length(esv) != N
        throw(DimensionMismatch(
            "the sectoral elasticity vector must have $N entries (got $(length(esv)))"))
    end
    sect = η0 || esv !== nothing
    if sect
        if init === nothing
            init = [ones(N); data.λ; ones(N); 0.0]
        elseif length(init) == 2N + 2
            init = [init[1:2N]; fill(init[2N+1], N); init[2N+2]]
        elseif length(init) == 2N + 1
            init = [init[1:2N]; fill(init[2N+1], N); 0.0]
        end
        length(init) == 3N + 1 || throw(DimensionMismatch(
            "the sectoral-wage closure expects a $(3N+1)-element init [p; y; w; F]"))
    elseif init === nothing
        # Default initialization: p=1, y=λ, wage=1, F=0
        init = [ones(N); data.λ; 1.0; 0.0]
    elseif length(init) == 2N + 1
        # Legacy init without the external transfer: F = 0
        init = [init; 0.0]
    end
    problem_f = sect ? problem_sectoral : problem
    # The sectoral system is far stiffer (Jacobian condition number ~9.4e7
    # against ~52.8 for the mobile all-N system, measured 2026-09-18): it starts
    # from a tighter Newton tolerance and gets a longer polish ladder, because
    # its external-account identity gate is 1e-12 (ADR-0020) and needs residuals
    # at ~1e-11 or below. The eta = 1 settings are unchanged (v5 behaviour).
    primary_tol = sect ? 1e-8 : 1e-6
    polish_steps = sect ? 6 : 4
    polish_tol = sect ? 1e-13 : 1e-10

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
            ProbN = NonlinearSolve.NonlinearProblem(problem_f, init, model)
            res = NonlinearSolve.solve(ProbN, reltol=primary_tol, abstol=primary_tol, maxiters=20000)
            # Quality gate = the ACTUAL residual, never the retcode (NonlinearSolve
            # reports Stalled on slow final convergence). Bounded LM polish (up to
            # 3 attempts) from the last point; verify before accepting. Mobile gate:
            # 1e-5 -- at the near-singular labour-equation direction the FD-Newton
            # floor is ~3e-6 (sum L off by 0.0003 percent, economically nil); the
            # budget identities remain EXACT and the DELTA equivalence gate stays
            # at machine precision.
            x = res.u
            rmax = maximum(abs, equilibrium_residuals(model, x))
            # ADR-0015: same monotone polish as `_solve_fixed` (see there), with
            # the eta = 0 target one decade tighter (`polish_tol` above).
            for _ in 1:polish_steps
                rmax <= polish_tol && break
                res = NonlinearSolve.solve(
                    NonlinearSolve.NonlinearProblem(problem_f, x, model),
                    NonlinearSolve.LevenbergMarquardt(); reltol=1e-12, abstol=1e-12, maxiters=20000)
                x_new = res.u
                r_new = maximum(abs, equilibrium_residuals(model, x_new))
                r_new < rmax || break
                x, rmax = x_new, r_new
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
    w = sect ? x[2N+1:3N] : x[2N+1]     # the sectoral wage vector
    F = sect ? x[3N+1] : x[2N+2]

    (; θ, ϵ, σ, η) = options.elasticities

    # Sectoral labor demand at equilibrium (the frozen allocation at eta = 0,
    # the cost-minimizing one at eta = 1)
    L_i = sectoral_labor_demand(p, q, w, model)

    # Wages vector (all equal to w at eta = 1; the sectoral vector at eta = 0)
    wages = sect ? collect(w) : fill(w, N)

    # Consumption — must match the budget-consistent, financed demand used
    # inside `problem`/`problem_sectoral` so the reported allocation is the
    # actual equilibrium allocation.
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    fin = model.financing
    ds_eff = preference_weights(fin, shocks.demand_shock)
    L_sum = sum(L_i)
    total_income = _wage_bill(w, L_i)
    E = household_expenditure(fin, model, total_income, p, L_sum; external_transfer = F)
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
    nominal_gdp = _wage_bill(w, L_i)

    return Solution(p, q, wages, consumption, numeraire, real_gdp_index, nominal_gdp, model; external_transfer = F)
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
function mobile_labor_model(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η::Float64; labor_bar::Union{Real, Nothing}=nothing, closure=:mobile, financing::Union{AbstractFinancing, Nothing}=nothing, eta_s::Union{Real, Nothing}=nothing, eta_s_vec::Union{AbstractVector, Nothing}=nothing)
    # BETA: an explicit supply elasticity forces the :beta closure. ADR-0022:
    # a sectoral elasticity vector selects the N-market system and likewise
    # belongs to :beta.
    if eta_s !== nothing || eta_s_vec !== nothing
        closure in (:mobile, :beta) || throw(ArgumentError(
            "eta_s / eta_s_vec apply to the :beta closure (got closure = $closure)"))
        closure = :beta
    end
    el = eta_s_vec !== nothing ?
        MobileLaborCESElasticities(θ, ϵ, σ, η, eta_s === nothing ? 0.0 : Float64(eta_s), eta_s_vec) :
        eta_s === nothing ? MobileLaborCESElasticities(θ, ϵ, σ, η) :
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
    labor_bar=nothing, closure=:mobile, financing=nothing, eta_s=nothing, eta_s_vec=nothing) = mobile_labor_model(data, shocks, Float64(θ), Float64(ϵ), Float64(σ), Float64(η);
    labor_bar=labor_bar === nothing ? nothing : Float64(labor_bar),
    closure=_closure_symbol(closure),
    financing=financing,
    eta_s=eta_s,
    eta_s_vec=eta_s_vec)
