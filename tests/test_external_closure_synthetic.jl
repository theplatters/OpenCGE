using OpenCGE, Test, LinearAlgebra

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the ADR-0019/ADR-0020 explicit external-account closure:
# every regime enforces all N goods-market clearings; the mobile eta = 1
# system carries the net external transfer F (entering household expenditure
# after tax); the eta = 0 endpoint is the SECTORAL-WAGE system (ADR-0020
# option C, a 3N+1 system), which closes the same identity by solving the
# per-sector FOC instead of pinning F; the canary
# `diff = S + T + M - (I+X) - (F + B_gov)` is the identity gap (approx 0 at
# every eta = 0/1 solution without legacy manna); F2/F3 are financing-neutral
# at both endpoints. Synthetic fixtures only: the real-table cells live with
# the paper repository. Tolerances via isapprox, never exact float equality.

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
		OpenCGE.problem(out_p, Xp, model)
		OpenCGE.problem(out_m, Xm, model)
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
