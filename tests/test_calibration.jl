using BeyondHulten
using DataFrames
using LinearAlgebra
using Test

# Contract tests for the Phase 3 calibration port (ADR-0006): DATASET_VARIANTS,
# drop_sectors, dataset_coverage, recalibrate_open, and the read_data datadir
# keyword. Fixture tests run everywhere; the real-data block is guarded on the
# (gitignored) IO table and skips when it is absent. Tolerances via isapprox,
# never exact float equality. The drop_sectors empty-components BoundsError
# below was a real bug: the first version of this test exposed it and
# src/core/calibration.jl now subsets value_added_components only when the row
# counts match.

"""Repo root (parent of tests/)."""
calibration_test_root() = normpath(joinpath(@__DIR__, ".."))

@testset "calibration: dataset variants" begin
    @test DATASET_VARIANTS["full"] == Int[]
    @test DATASET_VARIANTS["70s"] == [71]
    @test sort(DATASET_VARIANTS["reduced"]) == sort([71, 48, 18, 19, 53, 58, 13, 68])
end

@testset "calibration: drop_sectors on tiny_fixture" begin
    fx = tiny_fixture()
    # Empty drops are the identity (returns the object itself).
    @test drop_sectors(fx, Int[]) === fx
    d = drop_sectors(fx, [2])
    @test length(d.factor_share) == 1
    @test size(d.io, 1) == 1
    @test size(d.Ω) == (1, 1)
    @test size(d.Ω_raw) == (1, 1)
    # The v3 absorption vectors are subset with [keep] (ADR-0006 deviation),
    # not zeroed, so a dropped dataset keeps a valid household_baseline ...
    @test d.household_baseline == fx.household_baseline[[1]]
    @test d.gov_demand == fx.gov_demand[[1]]
    @test d.import_margin == fx.import_margin[[1]]
    @test d.exo_demand == fx.exo_demand[[1]]
    @test d.exports_demand == fx.exports_demand[[1]]
    @test d.factor_share == fx.factor_share[[1]]
    @test d.λ == fx.λ[[1]]
    # ... and the saving rate is kept.
    @test d.saving_rate == fx.saving_rate
end

@testset "calibration: dataset_coverage on tiny_fixture" begin
    fx = tiny_fixture()
    c = dataset_coverage(fx, [2])
    @test c.dropped_sectors == [2]
    @test c.gross_share_kept ≈ 50.0
    @test c.va_share_kept ≈ 50.0
    # fd_share_kept is NaN here by construction: the compact tiny fixture
    # carries zero domestic_final_demand (0/0), so it is not asserted.
    c0 = dataset_coverage(fx, Int[])
    @test c0.dropped_sectors == Int[]
    @test c0.gross_share_kept ≈ 100.0
    @test c0.va_share_kept ≈ 100.0
end

@testset "calibration: recalibrate_open argument validation" begin
    fx = tiny_fixture()
    @test_throws ArgumentError recalibrate_open(fx, "cbase2"; exo_scale = -0.1)
    @test_throws ArgumentError recalibrate_open(fx, "cbase2"; exo_scale = 1.5)
end

@testset "calibration: real-data v3 numbers (guarded)" begin
    root = calibration_test_root()
    io_path = joinpath(root, "data", "I-O_DE2019_formatiert.csv")
    if !isfile(io_path)
        @test_skip true  # data/ is gitignored; this block needs the IO table
    else
        full = read_data("I-O_DE2019_formatiert.csv"; datadir = root)
        @test length(full.factor_share) == 71
        recal71 = recalibrate_open(full, joinpath(root, "cbase2"); drops = Int[])
        @test recal71.saving_rate ≈ 0.3979 atol = 1e-4
        @test sum(recal71.gov_demand) ≈ 0.2141 atol = 1e-4
        dropped = drop_sectors(full, [71])
        @test length(dropped.factor_share) == 70
        recal70 = recalibrate_open(dropped, joinpath(root, "cbase2"); drops = [71])
        @test recal70.saving_rate ≈ 0.4259 atol = 1e-4
        @test sum(recal70.gov_demand) ≈ 0.2167 atol = 1e-4
        # Finiteness gate: worst-case round-gain column sums strictly below 1.
        colsums = (1.0 .- recal70.factor_share) .+
            (1.0 .- recal70.import_margin) .* (1.0 - recal70.saving_rate) .*
            recal70.factor_share
        @test maximum(colsums) < 1.0
        # Clamp mass recomputed from the returned fields (freeze record: 0.0873):
        # c0_dom = λ − Ω_raw'((1−fs)λ) − (1−m)(gG+inv) − expo, mass = Σ|min(c0,0)|.
        m = recal70.import_margin
        c0 = dropped.λ .- dropped.Ω_raw' * ((1.0 .- dropped.factor_share) .* dropped.λ) .-
            (1.0 .- m) .* (recal70.gov_demand .+ recal70.exo_demand) .- recal70.exports_demand
        @test sum(abs.(min.(c0, 0.0))) ≈ 0.0873 atol = 1e-4
    end
end
