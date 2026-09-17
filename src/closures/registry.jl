# src/closures/registry.jl — closure id → constructor mapping
#
# Small module-level mapping WITHOUT any status (ADR-0003). Status of each
# closure lives only in registry/closures.toml and renders on docs/status.md;
# this file only answers "which constructor implements id X?".
#
# Labour axis: BF / ALPHA / BETA / GAMMA / DELTA / ZETA (idea, no constructor).
# Financing axis: F1 / F2 / F3. Formulations: ADR-0002; registry: closures.toml.

const _CLOSURE_AXIS = Dict{Symbol,Symbol}(
    :BF => :labor, :ALPHA => :labor, :BETA => :labor,
    :GAMMA => :labor, :DELTA => :labor, :ZETA => :labor,
    :F1 => :financing, :F2 => :financing, :F3 => :financing,
)

"""All registered closure ids (labour + financing). `ZETA` maps to `nothing`."""
closure_ids() = [:BF, :ALPHA, :BETA, :GAMMA, :DELTA, :ZETA, :F1, :F2, :F3]

"""
	closure_axis(id) -> :labor | :financing

Axis of a registered closure id. Rejects unknown ids with `ArgumentError`.
"""
function closure_axis(id::Symbol)
    haskey(_CLOSURE_AXIS, id) || throw(ArgumentError(
        "unknown closure id $id; known ids are $(join(closure_ids(), ", "))"))
    _CLOSURE_AXIS[id]
end
closure_axis(id) = closure_axis(Symbol(id))

# ── Labour builders ──

"""
	bf_model(data, shocks, θ, ϵ, σ, η; kwargs...)

BF closure: geometric reallocation friction at `η`, flexible wage
(`closure=:mobile`). Extra kwargs (`labor_bar`, `financing`) pass through to
`mobile_labor_model`.
"""
bf_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real; kwargs...) =
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; closure = :mobile, kwargs...)

"""
	alpha_model(data, shocks, θ, ϵ, σ; kwargs...)

ALPHA closure: full-employment mobile labour — the η = 1 limit of BF by
construction (ADR-0002), flexible wage. Extra kwargs pass through.
"""
alpha_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real; kwargs...) =
    mobile_labor_model(data, shocks, θ, ϵ, σ, 1.0; closure = :mobile, kwargs...)

"""
	beta_model(data, shocks, θ, ϵ, σ, η; eta_s, kwargs...)

BETA closure: elastic total labour supply with elasticity `eta_s`
(`Σ L = L̄·(w/w0)^eta_s`). `eta_s` is required and forces the `:beta`
closure; `η` remains the BF reallocation parameter (never the supply
elasticity). Extra kwargs (`labor_bar`, `financing`) pass through.
"""
beta_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real;
        eta_s::Real, kwargs...) =
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; eta_s = eta_s, kwargs...)

"""
	gamma_model(data, shocks, θ, ϵ, σ, η; kwargs...)

GAMMA closure: fixed real wage (w = 1 numeraire), employment endogenous and
uncapped. Extra kwargs pass through (`labor_bar` must not be supplied).
"""
gamma_model(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real; kwargs...) =
    mobile_labor_model(data, shocks, θ, ϵ, σ, η; closure = :fixed, kwargs...)

"""
	closure_constructor(id) -> type/function/nothing

Constructor implementing a closure id: the labour builders above for
BF/ALPHA/BETA/GAMMA, `delta_model` for DELTA, the financing types for
F1/F2/F3, and `nothing` for ZETA (idea, no constructor). Rejects unknown
ids with `ArgumentError`.
"""
function closure_constructor(id::Symbol)
    id === :BF && return bf_model
    id === :ALPHA && return alpha_model
    id === :BETA && return beta_model
    id === :GAMMA && return gamma_model
    id === :DELTA && return delta_model
    id === :ZETA && return nothing
    id === :F1 && return PreferenceReallocation
    id === :F2 && return TaxFinanced
    id === :F3 && return ExternalDebt
    throw(ArgumentError(
        "unknown closure id $id; known ids are $(join(closure_ids(), ", "))"))
end
closure_constructor(id) = closure_constructor(Symbol(id))
