using BeyondHulten, Test, LinearAlgebra
using CSV, DataFrames

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for ADR-0022: the sectoral generalisation of BETA — N sectoral
# labour markets with a vector of real-wage supply elasticities, solved as the
# 3N+1 system `problem_sectoral`.
#
# The three claims under test:
#   1. the scalar BETA path is untouched (still the 2N+2 single-wage system);
#   2. the eta_s,i = 0 corner reproduces the ADR-0020 option C endpoint exactly
#      (the promotion's canary: the same equations, so the committed v6 BF cells
#      are the nested anchor);
#   3. a positive vector moves prices WITH the financing cell, while the
#      single-wage BETA closure at the same elasticity does not — the demand
#      channel is created by the N markets, not by the elasticity value.
# Real-table blocks skip when the (gitignored) IO or impulse table is absent.

"Repo root (parent of tests/)."
sec_test_root() = normpath(joinpath(@__DIR__, ".."))

"Full-71 A-bill calibration, or `nothing` when the IO table is absent."
function sec_test_data()
	root = sec_test_root()
	isfile(joinpath(root, "data", "I-O_DE2019_formatiert.csv")) || return nothing
	return recalibrate_open(
		read_data("I-O_DE2019_formatiert.csv"; datadir = root); exo_scale = 1.0)
end

"Matrix programme bundle `g` (and its incidence), or `nothing` when absent."
function sec_matrix_programme(data::Data)
	root = sec_test_root()
	imp_path = joinpath(root, "cbase2", "data_raw", "impulses.csv")
	isfile(imp_path) || return nothing
	imp = CSV.read(imp_path, DataFrame)
	rows = imp[imp.year .== 2024, :]
	raw = Matrix{Float64}(rows[1:1, 3:73])[:]
	ψ = raw ./ sum(raw)
	return (40300.0 / data.gdp_production) .* ψ, ψ
end

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
	@test_throws DimensionMismatch BeyondHulten.problem_sectoral(zeros(3N + 1),
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

@testset "ADR-0022 sectoral closure (full-71 calibration)" begin
	data = sec_test_data()
	if data === nothing
		@info "ADR-0022 real-data block skipped: IO table absent"
	else
		prog = sec_matrix_programme(data)
		if prog === nothing
			@info "ADR-0022 real-data block skipped: impulse table absent"
		else
			N = length(data.factor_share)
			shocks = Shocks(ones(N), ones(N), zeros(N))
			g, _ = prog

			# ── (2) THE CANARY: eta_s,i = 0 reproduces the eta = 0 endpoint ──
			# The committed matrix_5x3-v6 BF cells are exactly this system, so a
			# mismatch here means the kernel change is wrong.
			m0 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0;
				financing = TaxFinanced(g))
			s0 = solve(m0)
			mz = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
				eta_s_vec = zeros(N), financing = TaxFinanced(g))
			sz = solve(mz; init = [s0.prices_raw; s0.quantities; s0.wages_raw;
				s0.external_transfer])
			@test length(sz.wages_raw) == N
			@test maximum(abs, sz.prices_raw .- s0.prices_raw) < 1e-10
			@test maximum(abs, sz.quantities .- s0.quantities) < 1e-10
			@test maximum(abs, sz.wages_raw .- s0.wages_raw) < 1e-10
			@test abs(sz.external_transfer - s0.external_transfer) < 1e-10
			@test max_equilibrium_residual(sz) < 1e-9

			# ── (3) demand sensitivity, and the contrast with scalar BETA ──
			# Same eta_s value, two closures: the N-market system moves prices
			# with the demand composition; the single-wage one does not.
			eta = 0.5
			m_f2 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
				eta_s_vec = fill(eta, N), financing = TaxFinanced(g))
			m_f3 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
				eta_s_vec = fill(eta, N), financing = ExternalDebt(g))
			s_f2 = solve(m_f2; init = [s0.prices_raw; s0.quantities; s0.wages_raw;
				s0.external_transfer])
			s_f3 = solve(m_f3; init = [s_f2.prices_raw; s_f2.quantities; s_f2.wages_raw;
				s_f2.external_transfer])
			@test maximum(abs, s_f2.prices_raw .- 1) > 1e-3     # prices move at all
			@test maximum(abs, s_f2.prices_raw .- s_f3.prices_raw) < 1e-8  # F2/F3 neutral
			@test maximum(abs, s_f2.wages_raw .- s_f3.wages_raw) < 1e-8

			m_beta = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
				eta_s = eta, financing = TaxFinanced(g))
			s_beta = solve(m_beta; init = [s_f2.prices_raw; s_f2.quantities;
				s_f2.wages_raw[1]; s_f2.external_transfer])
			@test maximum(abs, s_beta.prices_raw .- 1) < 1e-8   # single wage: no price move

			# ── the external account closes in the sectoral system too ──
			X = [s_f2.prices_raw; s_f2.quantities; s_f2.wages_raw; s_f2.external_transfer]
			@test abs(external_balance_canary(m_f2, X).diff) < 1e-8
			@test max_equilibrium_residual(s_f2) < 1e-9
			@test gdp_components(m_f2, s_f2).wage_bill ≈
				sum(s_f2.wages_raw .* sectoral_labor_demand(s_f2.prices_raw, s_f2.quantities, s_f2.wages_raw, m_f2)) rtol = 1e-10
			@test gdp_wedge(s_f2) ≈ -external_balance_canary(m_f2, X).diff atol = 1e-9
		end
	end
end
