# src/closures/labor/types.jl — labor closure types (relocated from src/interface.jl / src/mobile_labor.jl)
#
# The canonical kernel defines exactly one set of closure descriptions.
# Status of each closure lives in registry/closures.toml (ADR-0003), never here.

"""Abstract supertype for labor-closure descriptions."""
abstract type AbstractLaborClosure end
"""Legacy CES closure, retaining its exogenous labor callback or symbol."""
struct ExogenousLaborClosure <: AbstractLaborClosure
	callback::Union{Function, Symbol}
end
"""Mobile labor with a flexible, market-clearing wage."""
struct FlexibleWageClosure <: AbstractLaborClosure end
"""Mobile labor with a fixed wage and unconstrained employment demand."""
struct FixedWageClosure <: AbstractLaborClosure end
"""
	SectoralElasticLaborClosure(eta_s_vec)

ADR-0022: the sectoral generalisation of BETA — N sectoral labour markets, each
with its own real-wage supply curve
`L^cm_i = Lbar_i * (w_i / CPI)^{eta_s,i}` (anchors `wbar_i / Pibar = 1`).
`eta_s_vec = zeros(N)` is exactly the eta = 0 sectoral-wage endpoint
(ADR-0020 option C); a uniform vector is BETA's supply rule applied market by
market. Selects the 3N+1 system `problem_sectoral`.
"""
struct SectoralElasticLaborClosure <: AbstractLaborClosure
	eta_s_vec::Vector{Float64}
end

"""
	ElasticLaborClosure(η_s; w0 = 1.0)

BETA closure: elastic total labour supply with elasticity `η_s` along the real
wage, anchored at the baseline wage `w0` (numeraire units). `η_s = 0` is
exactly ALPHA.
"""
struct ElasticLaborClosure <: AbstractLaborClosure
	η_s::Float64
	w0::Float64
	function ElasticLaborClosure(η_s::Real, w0::Real = 1.0)
		isfinite(η_s) && η_s ≥ 0 || throw(ArgumentError("η_s must be finite and non-negative"))
		w0 > 0 || throw(ArgumentError("the wage anchor w0 must be positive"))
		new(Float64(η_s), Float64(w0))
	end
end

CES(e::CESElasticities, c::ExogenousLaborClosure, reallocate::Bool=false) = CES(e, c.callback, reallocate)
labor_closure(options::CES) = ExogenousLaborClosure(options.labor_slack)
labor_closure(options::CobbDouglas) = ExogenousLaborClosure(options.labor_slack)

_closure_symbol(closure::Symbol) = closure
_closure_symbol(::FlexibleWageClosure) = :mobile
_closure_symbol(::FixedWageClosure) = :fixed
_closure_symbol(::ElasticLaborClosure) = :beta
_closure_symbol(::SectoralElasticLaborClosure) = :beta
_closure_symbol(closure) = throw(ArgumentError(
    "unsupported MobileLaborCES closure $closure; use :mobile, :fixed, :beta, FlexibleWageClosure(), FixedWageClosure(), or ElasticLaborClosure()"))
