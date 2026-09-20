using BeyondHulten, Test, LinearAlgebra
using CSV, DataFrames

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

"""Repo root (parent of tests/)."""
matrix_test_root() = normpath(joinpath(@__DIR__, ".."))

"""Two-sector fixture model/cell inputs shared by the headless testsets."""
function matrix_test_cells()
	data = tiny_fixture()
	model = mobile_labor_model(data, Shocks(ones(2), ones(2), zeros(2)),
		0.5, 0.5, 0.9, 1.0)
	ref = Solution([1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [0.5, 0.5],
		1.0, 1.0, 1.0, model)
	sol = Solution([1.1, 0.9], [1.2, 0.8], [1.0, 1.0], [0.6, 0.4],
		1.0, 1.0, 1.0, model)
	return data, model, ref, sol
end

@testset "matrix plots: canonical orders" begin
	@test matrix_labour_order() == ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
	@test matrix_financing_order() == ["F1", "F2", "F3"]
end

@testset "matrix plots: cell ids" begin
	ids = ["matrix_5x3-v10-BF-F2", "matrix_5x3-v10-BETA-F1-etas05",
		"matrix_5x3-v10-ALPHA-F1", "matrix_5x3-v10-GAMMA-F3",
		"matrix_5x3-v10-ZETA-F1"]
	@test matrix_cell_ids(ids, "matrix_5x3_v10") ==
		["matrix_5x3-v10-BF-F2", "matrix_5x3-v10-ALPHA-F1", "matrix_5x3-v10-GAMMA-F3"]
	# Canonical order, not input order: BF-F2 precedes ALPHA-F1.
	@test matrix_cell_ids(reverse(ids), "matrix_5x3_v10")[1] == "matrix_5x3-v10-BF-F2"
	@test matrix_cell_ids(String[], "matrix_5x3_v10") == String[]
	@test matrix_cell_ids(["other-design-BF-F1"], "matrix_5x3_v10") == String[]
end

@testset "matrix plots: baseline and cell data on the fixture" begin
	data, model, ref, sol = matrix_test_cells()
	base = matrix_baseline(data, ref)
	@test base.prices ≈ [1.0, 1.0]
	@test base.quantities ≈ [1.0, 1.0]
	@test base.wages ≈ [1.0, 1.0]
	@test base.consumption ≈ [0.5, 0.5]
	cell = matrix_cell_data("run-1", "ALPHA", "F2", 1.0, 0.5, "executed",
		data, sol, ref)
	@test cell.run_id == "run-1"
	@test cell.labour == "ALPHA"
	@test cell.financing == "F2"
	@test cell.eta == 1.0
	@test cell.eta_s == 0.5
	@test cell.status == "executed"
	@test cell.sectors == ["a", "b"]
	@test cell.labels == ["1", "2"]
	@test cell.prices ≈ [1.1, 0.9]
	@test cell.quantities ≈ [1.2, 0.8]
	@test cell.wages ≈ [1.0, 1.0]
	@test cell.consumption ≈ [1.2, 0.8]
	# Zero-baseline trade flows: zero level maps to 1.0.
	@test cell.imports ≈ [1.0, 1.0]
	@test cell.exports ≈ [1.0, 1.0]
	custom = matrix_cell_data("run-1", "ALPHA", "F2", 1.0, 0.5, "executed",
		data, sol, ref; labels = ["A", "B"],
		metrics = Dict("gdp" => 1.01), diagnostics = Dict("k" => 2))
	@test custom.labels == ["A", "B"]
	@test custom.metrics == Dict{String,Any}("gdp" => 1.01)
	@test custom.diagnostics == Dict{String,Any}("k" => 2)
	@test cell.metrics == Dict{String,Any}()
	# Identity cell: solution against its own baseline is all ones.
	self = matrix_cell_data("run-0", "BF", "F1", 0.0, 0.0, "executed",
		data, ref, ref)
	@test self.prices ≈ [1.0, 1.0]
	@test self.quantities ≈ [1.0, 1.0]
	@test self.wages ≈ [1.0, 1.0]
	@test self.consumption ≈ [1.0, 1.0]
end

@testset "matrix plots: dataset and frames" begin
	data, model, ref, sol = matrix_test_cells()
	base = matrix_baseline(data, ref)
	c1 = matrix_cell_data("id-BF-F2", "BF", "F2", 0.0, 0.0, "executed",
		data, sol, ref; metrics = Dict("gdp" => 1.01))
	c2 = matrix_cell_data("id-ALPHA-F1", "ALPHA", "F1", 1.0, 0.0, "executed",
		data, sol, ref; metrics = Dict("gdp" => 1.02, "employment" => 0.99))
	# Inserted out of canonical order on purpose.
	ds = matrix_dataset("matrix_5x3_v10", base, [c2, c1];
		provenance = Dict("src" => "test"), validation = Dict("id-BF-F2" => "ok"))
	@test ds.design == "matrix_5x3_v10"
	@test length(ds.cells) == 2
	@test ds.provenance == Dict{String,Any}("src" => "test")
	@test ds.validation == Dict("id-BF-F2" => "ok")
	frame = matrix_sectoral_frame(ds)
	@test names(frame) == ["run_id", "labour", "financing", "sector", "label", "variable", "rel"]
	@test nrow(frame) == 2 * 6 * 2
	# Canonical cell order (BF before ALPHA), then variable, then sector.
	@test frame.run_id[1] == "id-BF-F2"
	@test frame.run_id[end] == "id-ALPHA-F1"
	@test frame.variable[1:4] == ["price", "price", "quantity", "quantity"]
	@test frame.sector[1:2] == ["a", "b"]
	@test Set(frame.variable) ==
		Set(["price", "quantity", "wage", "consumption", "imports", "exports"])
	@test frame.rel[1] ≈ 1.1
	sumframe = matrix_summary_frame(ds)
	@test names(sumframe)[1:6] ==
		["run_id", "labour", "financing", "eta", "eta_s", "status"]
	@test "gdp" in names(sumframe) && "employment" in names(sumframe)
	@test nrow(sumframe) == 2
	@test sumframe.gdp ≈ [1.02, 1.01]
	@test sumframe.employment[1] ≈ 0.99
	@test isnan(sumframe.employment[2])
end

@testset "matrix plots: trade flows on the fixture" begin
	data, model, ref, sol = matrix_test_cells()
	flows = sectoral_trade_flows(model, sol)
	@test length(flows.imports) == 2
	@test length(flows.imports_final) == 2
	@test length(flows.imports_intermediate) == 2
	@test length(flows.exports) == 2
	@test flows.imports_final ≈ [0.0, 0.0] atol=1e-15
	@test flows.imports_intermediate ≈ [0.0, 0.0] atol=1e-15
	@test flows.exports ≈ [0.0, 0.0] atol=1e-15
	@test flows.imports ≈ flows.imports_final .+ flows.imports_intermediate
end

@testset "matrix plots: trade flows match gdp_components on the real table" begin
	root = matrix_test_root()
	io_path = joinpath(root, "data", "I-O_DE2019_formatiert.csv")
	imp_path = joinpath(root, "cbase2", "data_raw", "impulses.csv")
	if !isfile(io_path) || !isfile(imp_path)
		@test_skip true  # needs data/I-O_DE2019_formatiert.csv and cbase2 impulses.csv
	else
		data = recalibrate_open(read_data("I-O_DE2019_formatiert.csv"; datadir = root);
			exo_scale = 1.0)
		N = length(data.factor_share)
		imp = CSV.read(imp_path, DataFrame)
		rows = imp[imp.year .== 2024, :]
		raw = Matrix{Float64}(rows[1:1, 3:73])[:]
		ψ = raw ./ sum(raw)
		g = (40300.0 / data.gdp_production) .* ψ
		model = alpha_model(data, Shocks(ones(N), ones(N), zeros(N)),
			0.5, 0.5, 0.9; financing = TaxFinanced(g))
		sol = solve(model; init = [ones(N); data.λ; 1.0; 0.0])
		flows = sectoral_trade_flows(model, sol)
		comp = gdp_components(model, sol)
		@test isapprox(sum(flows.imports_final), -comp.V[5]; rtol = 1e-8, atol = 1e-10)
		@test isapprox(sum(flows.imports_intermediate), -comp.V[6]; rtol = 1e-8, atol = 1e-10)
		@test isapprox(sum(flows.exports), comp.V[4]; rtol = 1e-8, atol = 1e-10)
		@test flows.imports ≈ flows.imports_final .+ flows.imports_intermediate
	end
end

@testset "matrix plots: validate_cell" begin
	manifest = Dict("status" => "executed",
		"gates" => Dict("overall" => "pass"),
		"metrics" => Dict("gdp" => 1.01, "employment" => 0.99))
	re_metrics = Dict("gdp" => 1.01, "employment" => 0.99)
	rp = [1.0, 1.001]
	rq = [2.0, 1.999]
	ok = validate_cell(manifest, re_metrics, rp, rq, copy(rp), copy(rq))
	@test ok.pass
	@test isempty(ok.failures)
	@test ok.deltas["price"] ≈ 0.0 atol=1e-15
	@test ok.deltas["quantity"] ≈ 0.0 atol=1e-15
	@test ok.deltas["metric:gdp"] ≈ 0.0 atol=1e-15
	# Skipped stored arrays record NaN and still pass.
	skip = validate_cell(manifest, re_metrics, rp, rq, nothing, nothing)
	@test skip.pass
	@test isnan(skip.deltas["price"])
	@test isnan(skip.deltas["quantity"])
	# Failed status fails.
	bad_status = validate_cell(merge(manifest, Dict("status" => "failed")),
		re_metrics, rp, rq, nothing, nothing)
	@test !bad_status.pass
	@test any(occursin("status", f) for f in bad_status.failures)
	# Failed gate fails.
	bad_gate = validate_cell(merge(manifest, Dict("gates" => Dict("overall" => "fail"))),
		re_metrics, rp, rq, nothing, nothing)
	@test !bad_gate.pass
	@test any(occursin("overall", f) for f in bad_gate.failures)
	# Perturbed metric fails and names it.
	bad_metric = validate_cell(manifest, Dict("gdp" => 1.02, "employment" => 0.99),
		rp, rq, nothing, nothing)
	@test !bad_metric.pass
	@test any(f -> occursin("gdp", f), bad_metric.failures)
	@test bad_metric.deltas["metric:gdp"] ≈ 0.01
	# Stored price mismatch fails; wrong length fails.
	bad_price = validate_cell(manifest, re_metrics, rp, rq, [1.5, 1.001], copy(rq))
	@test !bad_price.pass
	@test any(f -> occursin("price", f), bad_price.failures)
	bad_len = validate_cell(manifest, re_metrics, rp, rq, [1.0], copy(rq))
	@test !bad_len.pass
	@test any(f -> occursin("length", f), bad_len.failures)
end

@testset "matrix plots: headless stubs require GLMakie" begin
	for f in (plot_matrix_overview, plot_matrix_wages, plot_matrix_prices,
			plot_matrix_quantities, plot_matrix_consumption, plot_matrix_trade,
			save_matrix_figures)
		err = try
			f()
			nothing
		catch e
			e
		end
		@test err isa ErrorException
		@test occursin("GLMakie", err.msg)
	end
end
