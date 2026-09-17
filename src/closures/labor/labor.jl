# ═══════════════════════════════════════════════════════════════════════════════
# Labor closures (Foundation I) — src/closures/labor/labor.jl
# ═══════════════════════════════════════════════════════════════════════════════
#
# Promoted from the frozen cbase2/src/closures.jl in Phase 2 (ADR-0005).
# Read cbase2 freely; never edit it. The design notes below are preserved
# from the cbase2 header as the docs pages (docs/closures/beta.md,
# docs/closures/delta.md); status lives only in registry/closures.toml
# (ADR-0003).
#
#   BETA  Elastic total labour supply (the referee-facing principal addition):
#         L^s = L̄ · [(w/P)/(w0/P0)]^{η_s},  interpreted as labour–leisure.
#         η_s = 0 reproduces ALPHA (vertical supply); η_s → ∞ approaches the
#         fixed-real-wage regime. Implemented as a labour-market-equation hook:
#         the :beta closure solves the SAME 2N+1 flexible-wage system as ALPHA
#         with the market-clearing equation Σ L_i = L̄ replaced by
#         Σ L_i = L̄ · (w/w0)^{η_s}. With the CPI numeraire (P = 1) the real
#         wage IS w; the anchor w0 = 1 (baseline numeraire wage).
#
#   DELTA IO endpoint (a corner of the matrix, never an independent row):
#         GAMMA (:fixed real factor price w = 1) + Leontief limit of the CES
#         core (θ, ϵ, σ → 0⁺, solved at a small ε for numerical robustness)
#         with full cost-minimizing allocation (η = 1). Under demand-only
#         shocks (A = 1) the zero-profit system pins prices at p = 1 exactly,
#         so relative prices are fixed and quantities absorb demand — the
#         Robinson (2006) fixed-price multiplier structure.
#
#   The analytic counterpart of the DELTA solve (exact linear Leontief system
#   with the endogenous consumption feedback) is provided by
#   `leontief_multiplier` and used as the equivalence test.
#
# `ElasticLaborClosure` itself lives in src/closures/labor/types.jl (with the
# other closure descriptions); the BETA residual hook and the DELTA/solve
# machinery live here. `ElasticLaborClosure` validation semantics are
# unchanged from cbase2 (η_s finite and non-negative, w0 positive).
# ═══════════════════════════════════════════════════════════════════════════════

"""
	labor_market_residual(closure, model, L_sum, w)

Total-labour-market residual of the flexible-wage systems. The BETA method
implements L^s = L̄ · (w/w0)^{η_s}; the ALPHA method lives in the core
(`src/core/equilibrium.jl`).
"""
function labor_market_residual(c::ElasticLaborClosure, model::Model{MobileLaborCES}, L_sum::Real, w::Real)
	L_sum - model.options.labor_bar * (max(w, eps(Float64)) / c.w0)^(c.η_s)
end

# ── DELTA corner ──────────────────────────────────────────────────────────────

"""
	delta_elasticities(ε = 1e-4)

Elasticities of the DELTA corner: the Leontief limit (θ, ϵ, σ → 0⁺, solved at
a small `ε` for numerical robustness) with full cost-minimizing allocation
(η = 1). Prices are demand-independent under `:fixed` + demand-only shocks, so
results are insensitive to `ε` up to O(ε).
"""
delta_elasticities(ε::Real = 1e-4) = MobileLaborCESElasticities(ε, ε, ε, 1.0)

"""
	delta_model(data, shocks; ε = 1e-4, financing = NoFinancing())

DELTA-corner model: `:fixed` closure + `delta_elasticities(ε)`.
"""
delta_model(data::Data, shocks::Shocks; ε::Real = 1e-4,
		financing::AbstractFinancing = NoFinancing()) =
	mobile_labor_model(data, shocks, delta_elasticities(ε).θ, delta_elasticities(ε).ϵ,
		delta_elasticities(ε).σ, delta_elasticities(ε).η; closure = :fixed, financing = financing)

"""
	leontief_multiplier(data, g; mode = :F3) -> (y, L, c_dom, F)

Analytic DELTA solution under the v2 open-absorption calibration: the exact
linear system of the model at (θ, ϵ, σ) = 0 under demand-only shocks
(A = 1 ⇒ p = 1, CPI = 1). With M = Ω_raw'·diag(1−fs), marginal consumption
gains (1−m)·ω per unit of (after-tax) wage income, employment L = fs'y, and
the reference income E_h0 = 1 − ΣgG:

  F3: y = M·y + (1−m)·ω·(1−τ0)·L + gG + (1−m)·g
  F2: y = M·y + (1−m)·ω·(L − ΣgG − Σg) + gG + (1−m)·g

The gain matrix has column sums strictly below 1 (import margins + the
proportional tax) — the v1 unit root is gone and the multiplier is finite.
`mode` selects the financing; the equivalence test asserts that
the nonlinear DELTA solve reproduces this system to O(ε).
"""
function leontief_multiplier(data::Data, g::Vector{Float64}; mode::Symbol = :F3)
	(; Ω_raw, factor_share, consumption_share, import_margin, gov_demand) = data
	n = length(factor_share)
	length(g) == n || throw(DimensionMismatch("g must have length $n"))
	mode in (:F2, :F3) || throw(ArgumentError("mode must be :F2 or :F3"))

	# Uniform-margin structure (matches the model): the household consumes
	# (1-s)E gross (saving sE leaks) with CPI-normalised weights omega
	# (consumption_share); only the domestic content (1 - m_i) circulates.
	# Exogenous injections: government gG (margin), investment I (margin),
	# exports X (NO margin). Income:
	#   E = (1 - tau0) * L            under F3 (programme untaxed; tau0 = sum gG)
	#   E = L - sum gG - sum g        under F2 (balanced-budget rule)
	M = Ω_raw' * Diagonal(1.0 .- factor_share)
	τ0 = sum(gov_demand)
	exo = (1.0 .- import_margin) .* (gov_demand + data.exo_demand) .+ data.exports_demand
	# F1/F3: FIXED NOMINAL baseline tax T = sum gG -> E = wL - sum gG, i.e. the
	# marginal after-tax share of employment income is 1 (the saving and import
	# leaks close the system). F2: balanced-budget rule E = L - sum gG - sum g.
	s̃ = 1.0
	gain = M + Diagonal(1.0 .- import_margin) * consumption_share * factor_share' * (1.0 - data.saving_rate) * s̃
	colsums = vec(sum(gain; dims=1))
	all(<(1), colsums) || error("finiteness gate failed: gain column sums must be < 1 (max = $(maximum(colsums)))")

	const_term = exo .+ (1.0 .- import_margin) .* g
	if mode === :F2
		const_term = const_term .- (1.0 .- import_margin) .* consumption_share .* (1.0 - data.saving_rate) .* (sum(gov_demand) + sum(g))
	else
		const_term = const_term .- (1.0 .- import_margin) .* consumption_share .* (1.0 - data.saving_rate) .* sum(gov_demand)
	end
	y = (I - gain) \ const_term   # equilibrium: y = G y + b  ⟺  (I − G) y = b
	L = dot(factor_share, y)
	# domestic household demand (for cross-check only; not used by the model)
	E = mode === :F3 ? (L - sum(gov_demand)) : (L - sum(gov_demand) - sum(g))
	c_dom = (1.0 .- import_margin) .* consumption_share .* (1.0 - data.saving_rate) .* E
	F = sum(import_margin .* g) + sum(import_margin .* (gov_demand + data.exo_demand)) +
		sum(import_margin .* consumption_share) * (1.0 - data.saving_rate) * E
	(; y = y, L = L, c_dom = c_dom, F = F)
end

"""
	solve_beta(data, shocks, θ, ϵ, σ, η; eta_s, financing = NoFinancing(),
	           init = nothing, steps = 5, max_refine = 8)

Standard solve strategy for the BETA closure: **η_s-continuation** from the
ALPHA equilibrium. The supply elasticity is raised along a ladder
`eta_s · k/steps`, warm-starting each step with the previous solution; the
last rung is `eta_s` itself. Direct single-shot solves at larger η_s can
require very many Newton iterations (the labour-demand/supply coupling is
stiff near the numeraire wage), while the homotopy converges in a handful of
cheap steps.

A rung is trusted only if its residuals verify below 1e-6; on failure the
rung is retried with (i) a Levenberg–Marquardt algorithm, (ii) a seeded
jitter of the warm start, and (iii) — failing both — a halved step toward
the previous rung (up to `max_refine` refinements). If refinement exhausts,
the continuation **jumps directly to the target** η_s from the last good
point: the BETA solution manifold can contain folds (recorded in
`cbase2/process_comments.md`), and a fold blocks step-wise progress while a
larger Newton step reaches the branch beyond it.
"""
function solve_beta(data::Data, shocks::Shocks, θ::Real, ϵ::Real, σ::Real, η::Real;
		eta_s::Real, financing::AbstractFinancing = NoFinancing(), init = nothing,
		steps::Int = 5, max_refine::Int = 8)
	steps ≥ 1 || throw(ArgumentError("steps must be ≥ 1"))
	target = Float64(eta_s)
	ladder = Float64[target * k / steps for k in 1:steps]
	x = init
	k = 1
	refine_count = 0
	while true
		η_k = ladder[k]
		mdl = mobile_labor_model(data, shocks, θ, ϵ, σ, η; financing = financing, eta_s = η_k)
		sol = _solve_rung_verified(mdl, x)
		if sol === nothing
			refine_count += 1
			if refine_count > max_refine
				# The refinement budget is exhausted: treat this as a fold of
				# the solution manifold and jump directly to the target.
				mdl_t = mobile_labor_model(data, shocks, θ, ϵ, σ, η;
					financing = financing, eta_s = target)
				sol_t = _solve_rung_verified(mdl_t, x)
				sol_t === nothing && error(
					"solve_beta: continuation stalled at η_s ≈ $(ladder[k-1]) and the jump to η_s = $target failed")
				return sol_t
			end
			prev = k == 1 ? 0.0 : ladder[k-1]
			insert!(ladder, k, (prev + η_k) / 2)
		else
			refine_count = 0
			x = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
			k == length(ladder) && return sol
			k += 1
		end
	end
end

# One continuation rung, solved with fallback strategies; returns nothing on
# failure. The QUALITY gate is the actual residual (< tol), never the solver
# retcode: NonlinearSolve reports `Stalled` when it stops making progress
# toward the requested tolerance even when the point already solves the
# system to well below the pipeline's 1e-6 residual gate — discarding such
# points was the cause of phantom continuation failures.
solve_verified(mdl::Model{MobileLaborCES}, x; tol::Float64 = 1e-6) =
	_solve_rung_verified(mdl, x; tol = tol)
function _solve_rung_verified(mdl::Model{MobileLaborCES}, x; tol::Float64 = 1e-6)
	for attempt in 1:4
		xx = try
			if attempt == 1
				# The core solve() settings (reltol 1e-6, maxiters 5000) are the
				# empirically proven configuration for this system — do not
				# "improve" them: tighter tolerances stall the solver.
				s = solve(mdl; init = x)
				[s.prices_raw; s.quantities; s.wages_raw[1]]
			elseif attempt == 2
				prob = NonlinearSolve.NonlinearProblem(problem, x, mdl)
				Float64.(NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt();
					reltol = 1e-6, abstol = 1e-6, maxiters = 1000).u)
			elseif attempt == 3
				# Seeded jitter: escape a stale convergence basin deterministically.
				jit = 1 .+ 1e-6 .* (cos.((1:length(x)) .* 1.7))
				prob = NonlinearSolve.NonlinearProblem(problem, x .* jit, mdl)
				Float64.(NonlinearSolve.solve(prob; reltol = 1e-6, abstol = 1e-6, maxiters = 1000).u)
			else
				prob = NonlinearSolve.NonlinearProblem(problem, x, mdl)
				Float64.(NonlinearSolve.solve(prob, NonlinearSolve.TrustRegion();
					reltol = 1e-6, abstol = 1e-6, maxiters = 1000).u)
			end
		catch
			continue
		end
		all(isfinite, xx) || continue
		maximum(abs, equilibrium_residuals(mdl, xx)) < tol || continue
		# Rebuild the Solution through the standard constructor path.
		return solve(mdl; init = xx)
	end
	return nothing
end
