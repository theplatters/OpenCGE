# ═══════════════════════════════════════════════════════════════════════════════
# Financing closures (Foundation II) — src/closures/financing/financing.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# Promoted from the frozen cbase2/src/financing.jl in Phase 2 (ADR-0005).
# Read cbase2 freely; never edit it. The design notes below are preserved
# from the cbase2 header as the docs pages (docs/closures/financing.md);
# status of F1/F2/F3 lives only in registry/closures.toml (ADR-0003).
#
# Design (agreed 2026-09-15, recorded in cbase2/process_comments.md):
#
#   Programme: real bundle g_i = m · G0 · ψ_i  (numeraire units), with ψ the
#   sectoral incidence (year-1 impulse shares) and G0 = 40,300 EUR m scaled
#   into model units by the notebook (1 model income unit = GDP at basic
#   prices ≈ 2,864,724 EUR m). Same ψ and scale for F1/F2/F3.
#
#   F1  PreferenceReallocation(d): compositional household demand shift.
#       Preference weights become demand_shock .* d; the household CES
#       normalizer keeps Σ p_i c_i = E exactly (budget closed by
#       construction). NO additive demand; at baseline prices the shift
#       redirects m·G0·ψ_i toward programme sectors and shrinks all other
#       categories pro rata within fixed E.
#
#   F2  TaxFinanced(g): additive real public/investment demand g_i;
#       endogenous public budget T(p) = Σ p_i g_i levied as a lump-sum tax:
#       household expenditure E = w·ΣL − T. Budget closes by construction
#       (Walras: Σ p(y − int − c − g) = wΣL − (wΣL − T) − T = 0).
#
#   F3  ExternalDebt(g): identical additive demand; household untaxed; the
#       external balance F = Σ p_i g_i is recorded post-solve (accounting
#       open-economy closure, import leakage mechanical via Ω_raw).
#
#   Compatibility (ADR-0005): the legacy unfinanced autonomous/investment
#   manna (Shocks.autonomous_demand / investment_shock) is RETAINED as an
#   explicitly documented compatibility path — cbase2 retired it, but the
#   root kernel keeps it so the characterization goldens and
#   rerun_results.jl/manuscript results stay reproducible. With NoFinancing
#   and zero v3 absorption fields the demand is numerically identical to the
#   pre-Phase-2 code. cbase2 experiments always pass zero manna vectors.
# ═══════════════════════════════════════════════════════════════════════════════

"""Baseline reference: no programme, household untaxed, no additive demand."""
struct NoFinancing <: AbstractFinancing end

"""
	PreferenceReallocation(shift)

F1 — compositional household demand shift. `shift` is a vector of strictly
positive preference shifters d_i (length N); the equilibrium budget closes
by construction through the household CES normalizer.
"""
struct PreferenceReallocation <: AbstractFinancing
	shift::Vector{Float64}
	function PreferenceReallocation(shift::Vector{<:Real})
		n = length(shift)
		n > 0 || throw(ArgumentError("preference shift vector must be non-empty"))
		all(>(0), shift) ||
			throw(ArgumentError("preference shifters must be strictly positive (a zero weight cannot be redirected)"))
		new(Float64.(shift))
	end
end

"""
	TaxFinanced(g)

F2 — additive real public/investment demand `g_i` (numeraire units),
financed by an endogenous lump-sum tax T(p) = Σ p_i g_i on household income.
"""
struct TaxFinanced <: AbstractFinancing
	g::Vector{Float64}
	function TaxFinanced(g::Vector{<:Real})
		n = length(g)
		n > 0 || throw(ArgumentError("public demand vector must be non-empty"))
		all(≥(0), g) ||
			throw(ArgumentError("public demand must be non-negative (use ExternalDebt/PreferenceReallocation for other designs)"))
		new(Float64.(g))
	end
end

"""
	ExternalDebt(g)

F3 — additive real public/investment demand `g_i` financed externally: the
household is untaxed and the external balance F = Σ p_i g_i records the
resource inflow post-solve.
"""
struct ExternalDebt <: AbstractFinancing
	g::Vector{Float64}
	function ExternalDebt(g::Vector{<:Real})
		n = length(g)
		n > 0 || throw(ArgumentError("public demand vector must be non-empty"))
		all(≥(0), g) ||
			throw(ArgumentError("public demand must be non-negative"))
		new(Float64.(g))
	end
end

# ── Demand-side hooks (called by problem() / problem_fixed() / solve()) ──

"""
	preference_weights(fin, demand_shock)

Effective household preference weights. F1 multiplies its shift into the
shock weights; all other closures pass the weights through unchanged.
"""
preference_weights(::AbstractFinancing, demand_shock::Vector{Float64}) = demand_shock
preference_weights(f::PreferenceReallocation, demand_shock::Vector{Float64}) =
	demand_shock .* f.shift

"""
	tau_rate(fin, model, p, w, L_sum)

Proportional income-tax rate (v3 open-economy calibration, homogeneous
budget): the government's real purchases gG are financed at CURRENT prices,
T = sum p_i gG_i, so the real tax burden is price-level-invariant and the
CPI numeraire selects the scale. F2 adds the programme bundle at the same
balanced rule: T = sum p_i (gG_i + g_i).
"""
tau_rate(::AbstractFinancing, model::Model{MobileLaborCES}, p, w::Real, L_sum::Real) =
	dot(p, model.data.gov_demand) / (w * L_sum)
tau_rate(::TaxFinanced, model::Model{MobileLaborCES}, p, w::Real, L_sum::Real) =
	dot(p, model.data.gov_demand .+ model.financing.g) / (w * L_sum)

"""
	household_expenditure(fin, model, wage_income, p, L_sum; external_transfer = 0.0)

Household expenditure base: after-tax wage income plus the net external
transfer, E = (1 - tau) w sum L + F. F is the equilibrium external margin of
the mobile system (ADR-0019): it enters AFTER tax (decision D1), so the
lump-sum tax base is wage income alone. The eta = 0 and fixed-wage systems
clear with F identically 0. The v2 import margin applies *within* spending
(the domestic content of each consumption category), not to E itself.
"""
household_expenditure(fin::AbstractFinancing, model::Model{MobileLaborCES},
		wage_income::Real, p, L_sum::Real; external_transfer::Real = 0.0) =
	(1 - tau_rate(fin, model, p, wage_income / max(L_sum, eps(Float64)), L_sum)) * wage_income + external_transfer

"""
	additive_demand(fin, N)

Gross programme bundle entering market clearing: g under F2/F3, zeros
otherwise (F1 shifts composition only). The DOMESTIC content applied by the
caller is `(1 .- data.import_margin) .* additive_demand(...)`.
"""
additive_demand(::AbstractFinancing, N::Int) = zeros(N)
additive_demand(f::Union{TaxFinanced, ExternalDebt}, N::Int) = f.g

"""
	public_budget(fin, p)

Government purchases to be financed domestically: baseline gG plus, under
F2, the programme bundle. Under F3 the programme is externally financed and
enters the external balance instead.
"""
public_budget(::AbstractFinancing, model::Model{MobileLaborCES}, p) =
	sum(model.data.gov_demand)
public_budget(f::TaxFinanced, model::Model{MobileLaborCES}, p) =
	sum(model.data.gov_demand) + dot(f.g, p)

"""
	external_balance(fin, model, p)

External balance: import content of the programme plus the import leakage of
the induced consumption response -- recorded post-solve under F3; zero
otherwise (F2's programme is domestically financed and the baseline
government is tax-financed).
"""
external_balance(::AbstractFinancing, model::Model{MobileLaborCES}, p) = 0.0
function external_balance(f::ExternalDebt, model::Model{MobileLaborCES}, p)
	dom_prog = (1 .- model.data.import_margin) .* f.g
	import_prog = dot(p, f.g) - dot(p, dom_prog)
	import_prog  # consumption-side leakage is reported separately by the notebooks
end

"""
	has_additive_anchor(fin)

Whether the financing carries an additive demand anchor that breaks the
homogeneity of the fixed-wage η ≈ 1 system (only F2/F3 do; F1 is purely
compositional and does NOT anchor scale — see the scale-indeterminacy
guard in `_solve_fixed`).
"""
has_additive_anchor(::AbstractFinancing) = false
has_additive_anchor(::Union{TaxFinanced, ExternalDebt}) = true
