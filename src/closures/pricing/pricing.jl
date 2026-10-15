# ═══════════════════════════════════════════════════════════════════════════════
# Pricing closures (Foundation III) — src/closures/pricing/pricing.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# ADR-0029 (the pricing door; full record in the paper repository). Demand can
# enter the price block only through the arguments of the unit-cost line,
#
#     cost = _ces_unit_cost(A_eff, factor_share, w, intermediate_price, ϵ),
#
# or as a markup on it. Three closures:
#
#   NoPricing()            the baseline: demand-invariant unit cost, p = c.
#   CapacityPricing(δ)     A_eff_i = A_i (y_i/λ_i)^(-δ): δ > 0 capacity
#                          pressure, δ < 0 the Kaldor–Verdoorn sign. Rent-free.
#   MarkupPricing(κ; s_pi) p_i = μ_i c_i, μ_i = 1 + κ (y_i/λ_i − 1). Creates a
#                          rent Π = Σ (1 − 1/μ_i) p_i y_i.
#
# `p` has two roles (ADR-0029 #4): the PURCHASER price that demand takes as
# given stays `p`; the multiplier of the first-order condition w = c ∂Φ/∂L is
# the MARGINAL cost `pFOC = p/μ`. Every FOC site reads `pFOC` and `A_eff`
# from the shared demand block (`_mobile_market_demand`), never
# `shocks.supply_shock` or `p` directly.
#
# The NoPricing hooks return their inputs unchanged (and `nothing` where the
# closure has no markup or rent), so the baseline is bit-identical by code
# structure rather than by multiplying with 1.0.
# ═══════════════════════════════════════════════════════════════════════════════

"""Baseline: demand-invariant unit cost, `p = c`, no rent."""
struct NoPricing <: AbstractPricing end

"""
	CapacityPricing(δ)

Route 1 (capacity / Kaldor–Verdoorn): effective technology
`A_eff_i = A_i (y_i/λ_i)^(-δ)`. `δ > 0` is capacity pressure (unit cost rises
with own output), `δ < 0` increasing returns. Rent-free; `δ = 0` nests the
baseline exactly.
"""
struct CapacityPricing <: AbstractPricing
	δ::Float64
	function CapacityPricing(δ::Real)
		isfinite(δ) || throw(ArgumentError("the capacity elasticity δ must be finite"))
		new(Float64(δ))
	end
end

"""
	MarkupPricing(κ; s_pi = 1.0)

Route 2 (distribution): `p_i = μ_i c_i` with `μ_i = 1 + κ (y_i/λ_i − 1)`. The
rent `Π = Σ (1 − 1/μ_i) p_i y_i` is retained (booked beside saving) at
`s_pi = 1` (Stage A); `s_pi < 1` distributes `(1 − s_pi) Π` to the household
(Stage B). `κ = 0` nests the baseline exactly.
"""
struct MarkupPricing <: AbstractPricing
	κ::Float64
	s_pi::Float64
	function MarkupPricing(κ::Real, s_pi::Real)
		isfinite(κ) || throw(ArgumentError("the markup slope κ must be finite"))
		0 <= s_pi <= 1 || throw(ArgumentError("the profit retention rate s_pi must lie in [0, 1]"))
		new(Float64(κ), Float64(s_pi))
	end
end
MarkupPricing(κ::Real; s_pi::Real = 1.0) = MarkupPricing(κ, s_pi)

# ── Hooks (one definition per expression) ──

"""
	pricing_technology(pr, A, y, λ) -> A_eff

Effective technology entering the unit cost and every FOC site.
"""
pricing_technology(::NoPricing, A, y, λ) = A

"""
	pricing_markup(pr, y, λ) -> μ

Markup vector, or `nothing` for a closure without a markup.
"""
pricing_markup(::NoPricing, y, λ) = nothing

"""
	pricing_foc_price(pr, p, μ) -> pFOC

Multiplier of the first-order conditions: the marginal cost `p/μ` (the
purchaser price `p` itself without a markup).
"""
pricing_foc_price(::NoPricing, p, μ) = p

"""
	pricing_profit(pr, p, μ, y) -> Π

Aggregate rent `Σ (1 − 1/μ_i) p_i y_i`, or `nothing` for a rent-free closure.
"""
pricing_profit(::NoPricing, p, μ, y) = nothing

"""
	pricing_bill_factor(pFOC, A_eff, ip, ϵ, ip_exponent = 1 - ϵ) -> k

The CES bill factor `k = pFOC^ϵ · A_eff^(ϵ−1) · ip^(ip_exponent)` (ADR-0016).
With the default exponent it values the intermediate-bill leaks `M_int` and
`T_int` (canary and `gdp_components`); with `θ − ϵ` it is the inner factor of
intermediate demand. The single definition of the FOC pattern `p^ϵ A^(ϵ−1)`.
"""
pricing_bill_factor(pFOC, A_eff, ip, ϵ, ip_exponent = 1 - ϵ) =
	pFOC .^ ϵ .* A_eff .^ (ϵ - 1) .* ip .^ ip_exponent

# Purchaser price from marginal cost: `μ .* c`, or `c` itself without a markup.
_apply_markup(::Nothing, c) = c

"""
	_pricing_state(model, p, y) -> (; A_eff, μ, pFOC)

The pricing hooks evaluated once at `(p, y)`; shared by `_mobile_market_demand`
and the public wrappers that are called outside it.
"""
function _pricing_state(model::Model{MobileLaborCES}, p, y)
	pr = model.pricing
	A_eff = pricing_technology(pr, model.shocks.supply_shock, y, model.data.λ)
	μ = pricing_markup(pr, y, model.data.λ)
	pFOC = pricing_foc_price(pr, p, μ)
	return (; A_eff, μ, pFOC)
end
