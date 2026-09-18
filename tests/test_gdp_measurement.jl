using BeyondHulten, Test, LinearAlgebra
using CSV, DataFrames

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the national-accounts GDP measurement API (ADR-0018):
# `gdp_components`, `gdp_deflator`, `gdp_income`, `gdp_expenditure`,
# `gdp_wedge` and `real_consumption` in src/core/diagnostics.jl. The
# real-table blocks need the (gitignored) IO table and the frozen impulse
# table; they skip when either file is absent. Tolerances via isapprox, never
# exact float equality.

"""Repo root (parent of tests/)."""
gdp_test_root() = normpath(joinpath(@__DIR__, ".."))

"""Full-71 A-bill calibration, or `nothing` when the IO table is absent."""
function gdp_test_data()
	root = gdp_test_root()
	isfile(joinpath(root, "data", "I-O_DE2019_formatiert.csv")) || return nothing
	return recalibrate_open(
		read_data("I-O_DE2019_formatiert.csv"; datadir = root); exo_scale = 1.0)
end

"""
Matrix programme bundle for `data`, or `nothing` when the impulse table is
absent. Mirrors the `programme_vectors` contract of experiments/run.jl for
the matrix_5x3_v3 `[programme]` section (total G0 = 40,300 EUR m scaled by
GDP at basic prices, renormalized 2024 impulse-share incidence): the ALPHA
cell anchors are incidence-independent at p = 1, but the GAMMA employment
response (0.0012598) is the impulse-incidence cell, not the uniform bundle.
"""
function gdp_matrix_programme(data::Data)
	root = gdp_test_root()
	imp_path = joinpath(root, "cbase2", "data_raw", "impulses.csv")
	isfile(imp_path) || return nothing
	imp = CSV.read(imp_path, DataFrame)
	rows = imp[imp.year .== 2024, :]
	raw = Matrix{Float64}(rows[1:1, 3:73])[:]
	ψ = raw ./ sum(raw)
	return (40300.0 / data.gdp_production) .* ψ
end

@testset "gdp measurement: real-table baseline identity" begin
	data = gdp_test_data()
	if data === nothing
		@test_skip true  # needs data/I-O_DE2019_formatiert.csv
	else
		N = length(data.factor_share)
		model = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)),
			0.5, 0.5, 0.9, 1.0)
		# The baseline is an exact root; solve returns immediately.
		sol = solve(model)
		c = gdp_components(model, sol)
		@test abs(c.wedge) < 1e-12
		@test c.wage_bill ≈ 1 atol=1e-12
		@test c.V ≈ [0.691675, 0.214101, 0.162367, 0.421945,
			-0.242941, -0.221403, -0.025745] atol=1e-6
		# ADR-0018 §3: all seven quantities are strictly positive here.
		@test all(>(0), c.Q)
	end
end

@testset "gdp measurement: known cell anchors" begin
	data = gdp_test_data()
	g = data === nothing ? nothing : gdp_matrix_programme(data)
	if data === nothing || g === nothing
		@test_skip true  # needs the IO table and cbase2 impulses.csv
	else
		N = length(data.factor_share)
		sh0 = Shocks(ones(N), ones(N), zeros(N))
		base = solve(mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0);
			init = [ones(N); data.λ; 1.0])
		# ALPHA-F2: demand-only programme, prices pinned at one.
		mA = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0;
			financing = TaxFinanced(g))
		solA = solve(mA; init = [ones(N); data.λ; 1.0])
		@test isapprox(real_consumption(solA) / real_consumption(base) - 1,
			-0.0169359; rtol = 1e-5)
		@test gdp_income(solA, base) - 1 ≈ 0 atol=1e-9
		@test gdp_deflator(solA, base) ≈ 1 atol=1e-9
		XA = [solA.prices_raw; solA.quantities; solA.wages_raw[1]]
		@test gdp_wedge(solA) ≈ -external_balance_canary(mA, XA).diff atol=1e-9
		# Fixed-wage GAMMA-F2 against the same mobile reference.
		mG = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0;
			closure = :fixed, financing = TaxFinanced(g))
		solG = solve(mG; init = [ones(N); data.λ])
		@test isapprox(gdp_income(solG, base) - 1, 0.0012598; rtol = 1e-5)
	end
end

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
	X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
	@test cs.wedge ≈ -external_balance_canary(m, X).diff atol=1e-9
end

@testset "gdp measurement: supply-shock deflator moves" begin
	data = gdp_test_data()
	if data === nothing
		@test_skip true  # needs data/I-O_DE2019_formatiert.csv
	else
		N = length(data.factor_share)
		sh0 = Shocks(ones(N), ones(N), zeros(N))
		base = solve(mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0);
			init = [ones(N); data.λ; 1.0])
		shS = Shocks(1 .+ 0.2 .* (1:N .== 1), ones(N), zeros(N))
		mS = mobile_labor_model(data, shS, 0.5, 0.5, 0.9, 1.0)
		solS = solve(mS; init = [ones(N); data.λ; 1.0])
		@test abs(gdp_deflator(solS, base) - 1) > 1e-6
		@test isfinite(gdp_income(solS, base)) && gdp_income(solS, base) > 0
		@test isfinite(gdp_expenditure(solS, base)) &&
			gdp_expenditure(solS, base) > 0
		XS = [solS.prices_raw; solS.quantities; solS.wages_raw[1]]
		@test gdp_wedge(solS) ≈ -external_balance_canary(mS, XS).diff atol=1e-5
	end
end
