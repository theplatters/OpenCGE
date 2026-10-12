using OpenCGE, Test, LinearAlgebra

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for ADR-0022: the sectoral generalisation of BETA - N sectoral
# labour markets with a vector of real-wage supply elasticities, solved as the
# 3N+1 system `problem_sectoral`. Structure-level claims only (synthetic
# fixtures): the scalar BETA path is untouched (still the 2N+2 single-wage
# system); the eta_s,i = 0 corner reproduces the ADR-0020 option C endpoint
# structurally; a wrong-length vector is rejected loudly; gdp_components uses
# the sectoral wage vector. The full-71 calibration block lives with the paper
# repository.

@testset "ADR-0022 sectoral labour markets (structure)" begin
	data = tiny_fixture()
	N = 2
	shocks = Shocks(ones(N), ones(N), zeros(N))

	# ── the scalar BETA path is untouched ──
	m_scalar = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; eta_s = 0.5)
	@test m_scalar.options.elasticities.eta_s_vec === nothing
	@test labor_closure(m_scalar.options) isa ElasticLaborClosure
	@test length(equilibrium_residuals(m_scalar, [ones(N); data.λ; 1.0; 0.0])) == 2N + 2

	# ── the sectoral path: 3N+1 unknowns, its own closure description ──
	m_sec = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; eta_s_vec = fill(0.5, N))
	@test m_sec.options.elasticities.eta_s_vec == fill(0.5, N)
	@test labor_closure(m_sec.options) isa SectoralElasticLaborClosure
	@test length(equilibrium_residuals(m_sec, [ones(N); data.λ; ones(N); 0.0])) == 3N + 1

	# ── a wrong-length vector is rejected loudly (never silently truncated) ──
	@test_throws DimensionMismatch solve(mobile_labor_model(
		data, shocks, 0.5, 0.5, 0.9, 1.0; eta_s_vec = fill(0.5, N + 1)))
	@test_throws DimensionMismatch OpenCGE.problem_sectoral(zeros(3N + 1),
		[ones(N); data.λ; ones(N); 0.0],
		mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; eta_s_vec = fill(0.5, N + 1)))

	# ── a zero vector is the eta = 0 endpoint: same closure, same length ──
	m_zero = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; eta_s_vec = zeros(N))
	@test length(equilibrium_residuals(m_zero, [ones(N); data.λ; ones(N); 0.0])) == 3N + 1

	# ── eta_s / eta_s_vec belong to the :beta closure only ──
	@test_throws ArgumentError mobile_labor_model(
		data, shocks, 0.5, 0.5, 0.9, 1.0; closure = :fixed, eta_s_vec = fill(0.5, N))

	# ── C1 regression: gdp_components must use the sectoral wage vector ──
	# The closed tiny_fixture makes the canary gap wage-independent, so the
	# open fixture (positive import margin + saving) is required to catch the
	# scalar-wage collapse on ADR-0022 cells.
	d0 = tiny_fixture()
	vals = Any[getfield(d0, f) for f in fieldnames(Data)]
	vals[findfirst(==(:import_margin), fieldnames(Data))] = [0.2, 0.2]
	vals[findfirst(==(:saving_rate), fieldnames(Data))] = 0.1
	d_open = Data(vals...)
	m_open = mobile_labor_model(d_open, Shocks(ones(N), ones(N), zeros(N)),
		0.5, 0.5, 0.9, 1.0; eta_s_vec = fill(0.5, N),
		financing = TaxFinanced([0.007, 0.003]))
	p = [1.10, 1.04]
	q = [0.48, 0.52]
	w = [1.20, 0.90]
	F = -0.02
	sol = Solution(p, q, w, [0.49, 0.51], 1.0, 1.0, 1.0, m_open; external_transfer = F)
	@test gdp_components(m_open, sol).wage_bill ≈ sum(w .* sectoral_labor_demand(p, q, w, m_open)) atol = 1e-12
	@test gdp_wedge(sol) ≈ -external_balance_canary(m_open, [p; q; w; F]).diff atol = 1e-12
end

