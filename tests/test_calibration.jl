using BeyondHulten
using DataFrames
using LinearAlgebra
using Test

# Contract tests for the Phase 3 calibration port (ADR-0006, ADR-0010):
# DATASET_VARIANTS, retained_io_table/retained_dataset (the review findings
# 2.2/2.3 repair of the old drop_sectors slicer), dataset_coverage,
# recalibrate_open, and the read_data datadir keyword. Fixture tests run
# everywhere; the real-data block is guarded on the (gitignored) IO table and
# skips when it is absent. Tolerances via isapprox, never exact float equality.

"""Repo root (parent of tests/)."""
calibration_test_root() = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

"""
3-sector IO table with the full Destatis row/column vocabulary, constructed so
that every `generate_data` accounting assertion holds and the retained slices
stay within the 10% GDP-reconciliation guard. Columns: Sektoren, the three
sector columns, the seven final-demand categories; rows: three sector rows,
the domestic-use row (73: "Gesamte Verwendung der inländischen Produktion",
the A-bill contract of ADR-0012), imports, goods taxes, and the five
value-added/production rows.
"""
function retained_io_fixture()
	sectors = ["A", "B", "C"]
	fdnames = [
		"Konsumausgaben der privaten Haushalte im Inland",
		"Konsumausgaben der privaten Organisationen o.E.",
		"Konsumausgaben des Staates",
		"Anlageinvestitionen f.Ausrüstungen u.sonst.Anlagen",
		"Anlageinvestitionen für Bauten",
		"Vorratsveränderungen und Nettozugang an Wertsachen",
		"Exporte",
	]
	Z = [10.0 5.0 2.0; 4.0 12.0 3.0; 2.0 3.0 8.0]   # Z[supplier, user]
	FD = [12.0 1.0 2.0 1.0 0.5 0.2 4.0;
	      15.0 1.0 2.0 0.5 0.3 0.1 5.0;
	       8.0 0.5 1.0 0.3 0.2 0.1 3.0]
	z = zeros(7)
	M = Matrix{Float64}(undef, 12, 10)
	for s in 1:3
		M[s, :] = vcat(Z[s, :], FD[s, :])
	end
	M[4, :]  = vcat([16.0, 20.0, 13.0], z)   # row 73: domestic intermediate bill (col sums of Z)
	M[5, :]  = vcat([2.0, 3.0, 1.0], [14.0, 0.9, 1.8, 1.1, 0.6, 0.2, 0.0])
	M[6, :]  = vcat([0.5, 0.5, 0.5], [9.0, 0.5, 1.2, 0.4, 0.2, 0.1, 0.0])
	M[7, :]  = vcat([10.0, 13.0, 6.5], z)
	M[8, :]  = vcat([6.0, 8.0, 4.0], z)
	M[9, :]  = vcat([1.0, 1.0, 0.5], z)
	M[10, :] = vcat([1.0, 1.5, 0.5], z)
	M[11, :] = vcat([2.0, 2.5, 1.5], z)
	M[12, :] = vcat([26.0, 33.0, 19.5], z)
	labels = vcat(sectors, ["Gesamte Verwendung der inländischen Produktion",
		"Verwendung der Importe",
		"Gütersteuern abzüglich Gütersubventionen", "Bruttowertschöpfung",
		"Arbeitnehmerentgelt im Inland",
		"Sonst.Produktionsabgaben abzgl. sonst.Subventionen",
		"Abschreibungen", "Nettobetriebsüberschuss", "Produktionswert"])
	df = DataFrame("Sektoren" => labels)
	for (j, name) in enumerate(vcat(sectors, fdnames))
		df[!, name] = M[:, j]
	end
	return df
end

"""Full-table `Data` built from the synthetic fixture table."""
function retained_fixture_data()
	io = retained_io_fixture()
	return BeyondHulten.assemble_data(io,
		BeyondHulten.generate_data(io; number_sectors = 3))
end

@testset "calibration: dataset variants" begin
	@test DATASET_VARIANTS["full"] == Int[]
	@test DATASET_VARIANTS["70s"] == [71]
	@test sort(DATASET_VARIANTS["reduced"]) == sort([71, 48, 18, 19, 53, 58, 13, 68])
end

@testset "calibration: final-demand columns located by name" begin
	io = retained_io_fixture()
	@test BeyondHulten.final_demand_columns(io) == collect(5:11)
	split = BeyondHulten.final_demand_split(io, 3)
	@test size(split.tot) == (3, 7)
	@test all(vec(sum(split.dom; dims = 1)) .<= vec(sum(split.tot; dims = 1)) .+ 1e-12)
end

@testset "calibration: retained_dataset on synthetic IO table" begin
	full = retained_fixture_data()
	# Empty drops are the identity (returns the object itself).
	@test retained_dataset(full, Int[]) === full

	d = retained_dataset(full, [2])
	@test length(d.factor_share) == 2
	@test size(d.Ω_raw) == (2, 2)
	@test size(d.io) == (11, 10)   # rows: 3 sectors + 9 aggregates − 1; cols: label + 3 sectors − 1 + 7 FD
	# Review findings 2.2/2.3: probability rows and the income unit survive.
	@test vec(sum(d.Ω_raw; dims = 2)) ≈ ones(2) atol=1e-12
	@test sum(d.labor_share) ≈ 1.0 atol=1e-12
	@test d.labor_share ≈ d.λ .* d.factor_share atol=1e-12
	# The legacy Törnqvist default is kept for an uncalibrated rebuild (ADR-0005).
	@test d.household_baseline ≈ d.consumption_share .* sum(d.labor_share) atol=1e-12

	# Recalibration runs on the rebuilt table and returns a normalized CPI block.
	de = recalibrate_open(d; exo_scale = 1.0)
	@test -1 < de.saving_rate < 1
	@test all(0 .<= de.import_margin .<= 1)
	@test sum(de.consumption_share) ≈ 1.0 atol=1e-12
	@test de.household_baseline ./ sum(de.household_baseline) ≈ de.consumption_share atol=1e-12
	@test all(>=(0), de.gov_demand) && all(>=(0), de.exo_demand) && all(>=(0), de.exports_demand)
	@test 1.0 - sum(de.gov_demand) > 0

	# Guards added with the repair (b46912d): loud failures instead of slicing.
	@test_throws ArgumentError retained_dataset(full, [0])
	@test_throws ArgumentError retained_dataset(full, [4])
	@test_throws ArgumentError retained_dataset(full, [1, 2, 3])
	io = retained_io_fixture()
	@test_throws ArgumentError retained_io_table(io, Int[]; number_sectors = 0)
	@test_throws ArgumentError retained_io_table(io, [3]; number_sectors = 2)
	@test_throws ArgumentError retained_io_table(io, [1, 2, 3]; number_sectors = 3)
end

@testset "calibration: dataset_coverage on synthetic fixture" begin
	full = retained_fixture_data()
	c = dataset_coverage(full, [2])
	@test c.dropped_sectors == [2]
	@test 0 < c.gross_share_kept < 100
	@test 0 < c.va_share_kept < 100
	@test 0 < c.fd_share_kept < 100
	c0 = dataset_coverage(full, Int[])
	@test c0.dropped_sectors == Int[]
	@test c0.gross_share_kept ≈ 100.0
	@test c0.va_share_kept ≈ 100.0
end

@testset "calibration: recalibrate_open argument validation" begin
	fx = tiny_fixture()
	@test_throws ArgumentError recalibrate_open(fx; exo_scale = -0.1)
	@test_throws ArgumentError recalibrate_open(fx; exo_scale = 1.5)
end

@testset "calibration: real-data v3 numbers (guarded)" begin
	root = calibration_test_root()
	io_path = joinpath(root, "data", "I-O_DE2019_formatiert.csv")
	if !isfile(io_path)
		@test_skip true  # data/ is gitignored; this block needs the IO table
	else
		full = read_data("I-O_DE2019_formatiert.csv"; datadir = root)
		@test length(full.factor_share) == 71
		@test sum(full.labor_share) ≈ 1.0 atol=1e-10
		recal71 = recalibrate_open(full; exo_scale = 1.0)
		# ADR-0012 (A-bill): the identity-implied saving rate; the clamp is gone.
		@test recal71.saving_rate ≈ 0.1199 atol = 1e-3
		@test sum(recal71.gov_demand) ≈ 0.214101 atol = 1e-4
		@test sum(recal71.M_int) ≈ 0.2214 atol = 1e-3
		# ADR-0013: row 75 (product taxes on intermediate use) is the third
		# component of the purchaser-price bill; booking it closes the canary.
		@test sum(recal71.T_int) ≈ 0.0257 atol = 1e-3
		@test sum(recal71.A_bill) + sum(recal71.M_int) + sum(recal71.T_int) ≈
			sum(recal71.λ) - 1.0 atol = 1e-10
		@test all(>=(0), recal71.household_baseline)

		# ADR-0013 contract: at the baseline root (and at the solved root) the
		# omitted N-th market residual equals the external canary to machine
		# precision. Before the T_int term the identity was short by exactly
		# sum(T_int) = 2.5745e-2 (the "−2.6e-2 reconciliation gap").
		mdl71 = mobile_labor_model(recal71, Shocks(ones(71), ones(71), zeros(71)),
			0.5, 0.5, 0.9, 1.0)
		X71 = [ones(71); recal71.λ; 1.0]
		@test abs(dot(ones(71), market_clearing_residuals(mdl71, X71)) -
			external_balance_canary(mdl71, X71).diff) < 1e-12
		sol71 = solve(mdl71)
		X71s = [sol71.prices_raw; sol71.quantities; sol71.wages_raw[1]]
		@test abs(dot(sol71.prices_raw, market_clearing_residuals(mdl71, X71s)) -
			external_balance_canary(mdl71, X71s).diff) < 1e-12

		# Review findings 2.2/2.3: the rebuilt 70-sector dataset has probability
		# rows (the old slice left 0.9713) and Σ labor_share = 1 (was 0.9878).
		dropped = retained_dataset(full, [71])
		@test length(dropped.factor_share) == 70
		@test vec(sum(dropped.Ω_raw; dims = 2)) ≈ ones(70) atol=1e-10
		@test sum(dropped.labor_share) ≈ 1.0 atol=1e-10
		recal70 = recalibrate_open(dropped; exo_scale = 1.0)
		# ADR-0012: 70s re-anchored rate; one microscopic identity residual
		# (−2.4e-6, the retained-economy gap) is clamped and documented.
		@test recal70.saving_rate ≈ 0.1285 atol = 1e-3
		@test sum(recal70.gov_demand) ≈ 0.216741 atol = 1e-4
		@test 1.0 - sum(recal70.gov_demand) ≈ 0.783259 atol = 1e-4
		@test sum(recal70.M_int) ≈ 0.2237 atol = 1e-3
		@test sum(recal70.T_int) ≈ 0.0260 atol = 1e-3
		# The 70s canary identity holds to machine precision around the
		# documented microscopic retained-economy residual (−2.4e-6), which
		# appears identically on both sides (ADR-0013).
		mdl70 = mobile_labor_model(recal70, Shocks(ones(70), ones(70), zeros(70)),
			0.5, 0.5, 0.9, 1.0)
		X70 = [ones(70); recal70.λ; 1.0]
		mk70 = dot(ones(70), market_clearing_residuals(mdl70, X70))
		@test abs(mk70 - external_balance_canary(mdl70, X70).diff) < 1e-12
		@test abs(mk70) < 1e-5   # the disclosed microscopic clamp
		# Finiteness gate: worst-case round-gain column sums strictly below 1.
		colsums = (1.0 .- recal70.factor_share) .+
			(1.0 .- recal70.import_margin) .* (1.0 - recal70.saving_rate) .*
			recal70.factor_share
		@test maximum(colsums) < 1.0
		# Clamp mass (review finding 2.5) is disclosed, not reconciled: recompute
		# from the returned fields and pin the modelled number.
		m = recal70.import_margin
		c0 = dropped.λ .- dropped.Ω_raw' * ((1.0 .- dropped.factor_share) .* dropped.λ) .-
			(1.0 .- m) .* (recal70.gov_demand .+ recal70.exo_demand) .- recal70.exports_demand
		@test sum(abs.(min.(c0, 0.0))) ≈ 0.086197 atol = 1e-4
	end
end
