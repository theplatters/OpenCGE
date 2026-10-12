using OpenCGE, Test, LinearAlgebra

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the national-accounts GDP measurement API (ADR-0018):
# `gdp_components`, `gdp_deflator`, `gdp_income`, `gdp_expenditure`,
# `gdp_wedge` and `real_consumption` in src/core/diagnostics.jl. Synthetic
# fixtures only: the real-table measurement anchors live with the paper
# repository. Tolerances via isapprox, never exact float equality.

@testset "gdp measurement: smoke fixture robustness" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	sh0 = Shocks(ones(N), ones(N), zeros(N))
	g = [0.007, 0.003]
	base_m = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0)
	base = solve(base_m)
	m = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0;
		financing = TaxFinanced(g))
	sol = solve(m; init = [ones(N); data.λ; 1.0])
	cb = gdp_components(base_m, base)
	cs = gdp_components(m, sol)
	# The closed fixture carries no gov/investment/exports/intermediate
	# leaks, so those components are zero at base or current; the composite
	# must handle the zero-base/zero-current slots without dividing by zero.
	@test cb.V[2] == 0 && cs.V[2] != 0
	@test all(==(0), cb.V[3:4]) && all(==(0), cs.V[3:4])
	@test all(==(0), cb.V[6:7]) && all(==(0), cs.V[6:7])
	@test isfinite(gdp_deflator(sol, base)) && gdp_deflator(sol, base) > 0
	@test isfinite(gdp_income(sol, base)) && gdp_income(sol, base) > 0
	@test isfinite(gdp_expenditure(sol, base)) && gdp_expenditure(sol, base) > 0
	@test isfinite(gdp_wedge(sol))
	@test real_consumption(sol) == sol.real_gdp
	# The closed cores have no open-economy blocks.
	@test_throws ArgumentError gdp_components(Model(data, sh0, CES()), sol)
	# The wedge/canary identity holds on the fixture as well.
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
	@test cs.wedge ≈ -external_balance_canary(m, X).diff atol=1e-9
end

