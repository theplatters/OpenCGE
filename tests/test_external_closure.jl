using BeyondHulten, Test, LinearAlgebra
using CSV, DataFrames

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the ADR-0019/ADR-0020 explicit external-account closure:
# every regime enforces all N goods-market clearings; the mobile η = 1 system
# carries the net external transfer F (entering household expenditure after
# tax); the η = 0 endpoint is the SECTORAL-WAGE system (ADR-0020 option C, a
# 3N+1 system), which closes the same identity by solving the per-sector FOC
# instead of pinning F; the canary `diff = S + T + M − (I+X) − (F + B_gov)` is
# the identity gap (≈ 0 at every η = 0/1 solution without legacy manna); F2/F3
# are financing-neutral at both endpoints.
# Real-table blocks need the (gitignored) IO table and the frozen impulse
# table; they skip when either file is absent. Tolerances via isapprox, never
# exact float equality.

"""Repo root (parent of tests/)."""
ext_test_root() = normpath(joinpath(@__DIR__, ".."))

"""Full-71 A-bill calibration, or `nothing` when the IO table is absent."""
function ext_test_data()
	root = ext_test_root()
	isfile(joinpath(root, "data", "I-O_DE2019_formatiert.csv")) || return nothing
	return recalibrate_open(
		read_data("I-O_DE2019_formatiert.csv"; datadir = root); exo_scale = 1.0)
end

"""
Matrix programme bundle for `data`, or `nothing` when the impulse table is
absent. Same incidence as tests/test_gdp_measurement.jl (`gdp_matrix_programme`):
total G0 = 40,300 EUR m scaled by GDP at basic prices with the renormalized
2024 impulse-share incidence.
"""
function ext_matrix_programme(data::Data)
	root = ext_test_root()
	imp_path = joinpath(root, "cbase2", "data_raw", "impulses.csv")
	isfile(imp_path) || return nothing
	imp = CSV.read(imp_path, DataFrame)
	rows = imp[imp.year .== 2024, :]
	raw = Matrix{Float64}(rows[1:1, 3:73])[:]
	ψ = raw ./ sum(raw)
	return (40300.0 / data.gdp_production) .* ψ, ψ
end

"""
F1 preference-tilt weights. Replicates the documented `d = 1 .+ G0 .* ψ1 ./
max.(c0, 1e-12)` formula (ψ restricted to positive-baseline sectors and
renormalized) whose canonical implementation is `f1_tilt_weights` in
experiments/run.jl (not included here by design).
"""
function ext_tilt(c0::AbstractVector, ψ::AbstractVector, g::AbstractVector)
	pos = c0 .> 0
	ψ1 = (ψ .* pos) ./ sum(ψ .* pos)
	return 1.0 .+ sum(g) .* ψ1 ./ max.(c0, 1e-12)
end

"""Financing closures for the matrix cells (F1 tilt, F2 tax, F3 external debt)."""
function ext_financing(data::Data, g::Vector{Float64}, ψ::Vector{Float64})
	return (
		F1 = PreferenceReallocation(ext_tilt(data.household_baseline, ψ, g)),
		F2 = TaxFinanced(g),
		F3 = ExternalDebt(g),
	)
end

@testset "external closure: real-table baseline reproduction" begin
	data = ext_test_data()
	if data === nothing
		@test_skip true  # needs data/I-O_DE2019_formatiert.csv
	else
		N = length(data.factor_share)
		model = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)),
			0.5, 0.5, 0.9, 1.0)
		# The baseline is an exact root with F = 0; solve returns immediately.
		sol = solve(model)
		X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
		@test maximum(abs, equilibrium_residuals(model, X)) ≤ 1e-12
		@test abs(sol.external_transfer) ≤ 1e-12
		@test sol.prices_raw ≈ ones(N) atol=1e-10
		@test sol.quantities ≈ data.λ atol=1e-10
	end
end

@testset "external closure: real-table ALPHA cells" begin
	data = ext_test_data()
	pg = data === nothing ? nothing : ext_matrix_programme(data)
	if data === nothing || pg === nothing
		@test_skip true  # needs the IO table and cbase2 impulses.csv
	else
		N = length(data.factor_share)
		g, ψ = pg
		fins = ext_financing(data, g, ψ)
		sh0 = Shocks(ones(N), ones(N), zeros(N))
		base = solve(mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0);
			init = [ones(N); data.λ; 1.0])
		sols = Dict{Symbol,Any}()
		cans = Dict{Symbol,Any}()
		# Measured 2026-09-18: F = +7.848e-4 (F1) / −1.014e-3 (F2) /
		# −1.4324e-2 (F3); canary gap |diff| ≤ 3e-16 at all three.
		for (key, F_pin) in ((:F1, 7.848e-4), (:F2, -1.014e-3), (:F3, -1.4324e-2))
			m = alpha_model(data, sh0, 0.5, 0.5, 0.9; financing = fins[key])
			sol = solve(m; init = [ones(N); data.λ; 1.0])
			X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
			@test maximum(abs, equilibrium_residuals(m, X)) < 1e-10
			@test maximum(abs, market_clearing_residuals(m, X)) ≤ 1e-10
			can = external_balance_canary(m, X)
			@test abs(can.diff) ≤ 1e-12
			@test sol.external_transfer ≈ F_pin atol=1e-6
			# The booked external position is recorded, not lost.
			@test can.financing ≈ sol.external_transfer + can.programme_financing atol=1e-12
			@test gdp_income(sol, base) ≈ gdp_expenditure(sol, base) atol=1e-12
			sols[key] = sol
			cans[key] = can
		end
		@test cans[:F1].programme_financing ≈ 0 atol=1e-12
		@test cans[:F2].programme_financing ≈ 0 atol=1e-12
		@test cans[:F3].programme_financing ≈ 1.3310e-2 atol=1e-6
		# FINANCING NEUTRALITY (theorem of the closure): F2 and F3 have
		# identical real equilibria with F_F3 = F_F2 − B_gov, hence identical
		# net external positions F + B_gov.
		s2, s3 = sols[:F2], sols[:F3]
		@test maximum(abs, s2.prices_raw .- s3.prices_raw) ≤ 1e-10
		@test maximum(abs, s2.quantities .- s3.quantities) ≤ 1e-10
		@test abs(s2.wages_raw[1] - s3.wages_raw[1]) ≤ 1e-10
		@test abs(s3.external_transfer -
			(s2.external_transfer - cans[:F3].programme_financing)) ≤ 1e-10
		@test abs(cans[:F2].financing - cans[:F3].financing) ≤ 1e-12
	end
end

@testset "external closure: real-table BF cells (sectoral wages, closed account)" begin
	data = ext_test_data()
	pg = data === nothing ? nothing : ext_matrix_programme(data)
	if data === nothing || pg === nothing
		@test_skip true  # needs the IO table and cbase2 impulses.csv
	else
		N = length(data.factor_share)
		g, ψ = pg
		fins = ext_financing(data, g, ψ)
		sh0 = Shocks(ones(N), ones(N), zeros(N))
		# ADR-0020 option C: at η = 0 the endpoint solves the per-sector FOC, so
		# the frozen allocation is cost-minimizing at the sectoral wages, the
		# identity gap vanishes and F is identified. The retired F = 0 pin left F
		# unidentified (the reported position collapsed to B_gov = Σp·g).
		# Measured 2026-09-18 on the full-71 table: F = −5.815e-3 (F1) /
		# −8.528e-3 (F2) / −2.337e-2 (F3); booked = −5.815e-3 / −8.528e-3 /
		# −8.528e-3; |gap| ≤ 1.04e-11 (BF-F1, the stiffest cell).
		sols = Dict{Symbol,Any}()
		cans = Dict{Symbol,Any}()
		for (key, F_ref) in ((:F1, -0.0058149704), (:F2, -0.0085284060), (:F3, -0.0233746904))
			m = bf_model(data, sh0, 0.5, 0.5, 0.9, 0.0; financing = fins[key])
			sol = solve(m)
			@test length(sol.wages_raw) == N
			X = [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer]
			@test length(X) == 3N + 1
			@test maximum(abs, equilibrium_residuals(m, X)) < 1e-10
			@test maximum(abs, market_clearing_residuals(m, X)) ≤ 1e-10
			@test sectoral_labor_gap(m, sol.prices_raw, sol.quantities,
				sol.wages_raw) ≤ 1e-10
			can = external_balance_canary(m, X)
			@test abs(can.diff) ≤ 1e-10        # the account closes at η = 0
			@test sol.external_transfer ≈ F_ref atol=1e-6
			@test can.financing ≈ sol.external_transfer + can.programme_financing atol=1e-12
			sols[key] = sol
			cans[key] = can
		end
		# The sectoral wage vector is a real instrument, not a rescaling of one
		# common wage: measured range 0.9702 … 1.6719 on the full-71 table.
		@test minimum(sols[:F3].wages_raw) < 0.99
		@test maximum(sols[:F3].wages_raw) > 1.6
		@test cans[:F3].programme_financing ≈ 1.4846e-2 atol=1e-6
		# FINANCING NEUTRALITY at η = 0 (as at η = 1): identical real equilibria
		# with F_F3 = F_F2 − B_gov, hence identical net external positions.
		s2, s3 = sols[:F2], sols[:F3]
		@test maximum(abs, s2.prices_raw .- s3.prices_raw) ≤ 1e-9
		@test maximum(abs, s2.quantities .- s3.quantities) ≤ 1e-9
		@test maximum(abs, s2.wages_raw .- s3.wages_raw) ≤ 1e-9
		@test abs(s3.external_transfer -
			(s2.external_transfer - cans[:F3].programme_financing)) ≤ 1e-9
		@test abs(cans[:F2].financing - cans[:F3].financing) ≤ 1e-9
	end
end

@testset "external closure: fixed-wage GAMMA cell" begin
	data = ext_test_data()
	pg = data === nothing ? nothing : ext_matrix_programme(data)
	if data === nothing || pg === nothing
		@test_skip true  # needs the IO table and cbase2 impulses.csv
	else
		N = length(data.factor_share)
		g, ψ = pg
		sh0 = Shocks(ones(N), ones(N), zeros(N))
		fin = TaxFinanced(g)
		m = gamma_model(data, sh0, 0.5, 0.5, 0.9, 1.0; financing = fin)
		sol = solve(m)
		@test maximum(abs,
			equilibrium_residuals(m, [sol.prices_raw; sol.quantities])) ≤ 1e-6
		@test abs(sol.external_transfer) ≤ 1e-12
		can = external_balance_canary(m, [sol.prices_raw; sol.quantities])
		@test abs(can.diff) ≤ 1e-9
	end
end

@testset "external closure: closed fixture solves with F = 0" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
	sol = solve(model)
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
	@test maximum(abs, equilibrium_residuals(model, X)) ≤ 1e-12
	@test abs(sol.external_transfer) ≤ 1e-12
end

@testset "external closure: F is identified (non-degenerate column)" begin
	# Finite-difference Jacobian of `problem` at the closed-fixture solution:
	# the F column must carry real variation (F is an identified unknown, not
	# a null direction). Measured F-column norm 0.71 on 2026-09-18.
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
	sol = solve(model)
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
	n = length(X)
	h = 1e-7
	Fcol = let
		out_p = similar(X)
		out_m = similar(X)
		Xp = copy(X)
		Xp[end] += h
		Xm = copy(X)
		Xm[end] -= h
		BeyondHulten.problem(out_p, Xp, model)
		BeyondHulten.problem(out_m, Xm, model)
		(out_p .- out_m) ./ (2h)
	end
	@test norm(Fcol) > 1e-3
end

@testset "external closure: multi-start invariance on the fixture" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
	X0 = [ones(N); data.λ; 1.0; 0.0]
	jitter = 1e-4 * collect(1.0:length(X0))   # small deterministic perturbation
	s1 = solve(model; init = X0)
	s2 = solve(model; init = X0 .+ jitter)
	@test maximum(abs, s1.prices_raw .- s2.prices_raw) ≤ 1e-8
	@test maximum(abs, s1.quantities .- s2.quantities) ≤ 1e-8
	@test abs(s1.external_transfer - s2.external_transfer) ≤ 1e-8
end

@testset "external closure: eta = 0 sectoral wages on the fixture" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	model = bf_model(data, shocks, 0.5, 0.5, 0.9, 0.0)
	sol = solve(model)
	X = [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer]
	@test length(X) == 3N + 1
	# The closed fixture carries no external block, so the identity forces F = 0
	# and the sectoral FOC is met: the account closes exactly. This fixture's
	# baseline is self-consistent (q = λ), so the sectoral wages stay at 1.
	@test abs(sol.external_transfer) ≤ 1e-12
	@test maximum(abs, equilibrium_residuals(model, X)) ≤ 1e-10
	@test maximum(abs, market_clearing_residuals(model, X)) ≤ 1e-10
	@test sectoral_labor_gap(model, sol.prices_raw, sol.quantities,
		sol.wages_raw) ≤ 1e-10
	@test abs(external_balance_canary(model, X).diff) ≤ 1e-10
	@test sol.wages_raw ≈ ones(N) atol=1e-6
end

@testset "external closure: legacy 2N+1 vectors append F = 0" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
	sol = solve(model)
	p, y, w = sol.prices_raw, sol.quantities, sol.wages_raw[1]
	@test equilibrium_residuals(model, [p; y; w]) ≈
		equilibrium_residuals(model, [p; y; w; 0.0])
end

@testset "external closure: exact book identity on the fixture" begin
	# `gdp_components(...).wedge ≈ −canary.diff` on and off equilibrium, and
	# `external_financing ≈ external_transfer + programme_financing`. The F3
	# cell exercises a nonzero programme booking (measured B_gov = +0.01 with
	# F = −0.01 on the closed fixture).
	data = tiny_fixture()
	N = length(data.factor_share)
	shocks = Shocks(ones(N), ones(N), zeros(N))
	g = [0.007, 0.003]
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0;
		financing = ExternalDebt(g))
	sol = solve(model)
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
	c = gdp_components(model, sol)
	can = external_balance_canary(model, X)
	@test c.wedge ≈ -can.diff atol=1e-9
	@test c.external_financing ≈ c.external_transfer + c.programme_financing atol=1e-12
	@test can.programme_financing ≈ 0.01 atol=1e-12
	# Off equilibrium: the same identity at an arbitrary point with F = 0.01.
	off = Solution(ones(N), data.λ, ones(N), sol.consumption, 1.0, 1.0, 1.0,
		model; external_transfer = 0.01)
	Xoff = [ones(N); data.λ; 1.0; 0.01]
	coff = gdp_components(model, off)
	canoff = external_balance_canary(model, Xoff)
	@test coff.wedge ≈ -canoff.diff atol=1e-9
	@test coff.external_financing ≈
		coff.external_transfer + coff.programme_financing atol=1e-12
end

@testset "external closure: legacy manna is absorbed by F and reported" begin
	# ADR-0005 compatibility path: nonzero autonomous/investment manna is
	# unfinanced demand that is not booked in B_gov. At an all-N η = 1
	# solution the identity gap reads exactly the manna's value p·(A+G) and
	# the free F absorbs it (F = −p·(A+G) at zero programme); all clearings
	# still hold, and the wedge duality stands. The matrix designs pass zero
	# manna, where the gap is ≈ 0 (see the other testsets).
	data = tiny_fixture()
	N = length(data.factor_share)
	autonomous = [0.1, 0.0]
	shocks = Shocks(ones(N), ones(N); autonomous_demand = autonomous)
	model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
	sol = solve(model)
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
	can = external_balance_canary(model, X)
	A = autonomous .* data.consumption_share .* sum(data.labor_share)
	manna_val = dot(sol.prices_raw, A)
	@test maximum(abs, market_clearing_residuals(model, X)) ≤ 1e-10
	@test manna_val > 1e-6                    # the fixture carries manna
	@test can.diff ≈ manna_val atol=1e-9      # the gap is the unbooked manna
	@test sol.external_transfer ≈ -manna_val atol=1e-9
	@test gdp_components(model, sol).wedge ≈ -can.diff atol=1e-9
end
