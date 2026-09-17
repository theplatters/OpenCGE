using BeyondHulten
using Test

# Contract tests for the F1 preference tilt (ports cbase2/scripts/verify_v3.jl
# lines ~146-160, revisefinal commit 0f33ad6; cbase2/review.md defect 9).
# The retired `1 .+ ψ` stand-in was a non-reference tilt (ψ added directly
# instead of `G0 .* ψ1 ./ c0`); the
# canonical tilt is `d = 1 .+ G0 .* ψ1 ./ c0` with ψ restricted to
# positive-baseline sectors and renormalised.

# run.jl is include-safe (PROGRAM_FILE guard); skip the re-include when another
# test file already loaded it (both are included into Main by runtests.jl).
isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))
isdefined(Main, :run_design) || include(joinpath(@__DIR__, "..", "experiments", "run.jl"))

@testset "f1 tilt: formula against a hand computation" begin
    data = tiny_fixture()
    @test data.household_baseline ≈ [0.5, 0.5]
    ψ = [0.7, 0.3]
    g = [0.07, 0.03]
    d = f1_tilt_weights(data.household_baseline, ψ, g)
    G0 = sum(g)
    expected = 1.0 .+ G0 .* ψ ./ data.household_baseline
    @test d ≈ expected
    # Hand values: G0 = 0.1, d = [1.14, 1.06].
    @test d ≈ [1.14, 1.06]
end

@testset "f1 tilt: ψ restricted to positive-baseline sectors and renormalised" begin
    baseline = [0.4, 0.0, -0.1, 0.6]
    ψ = [0.25, 0.25, 0.25, 0.25]
    g = [0.025, 0.025, 0.025, 0.025]
    d = f1_tilt_weights(baseline, ψ, g)
    G0 = sum(g)
    # Only sectors 1 and 4 carry mass: ψ1 = [0.5, 0, 0, 0.5].
    @test d[1] ≈ 1.0 + G0 * 0.5 / 0.4
    @test d[4] ≈ 1.0 + G0 * 0.5 / 0.6
    # Zero/negative-baseline sectors are untouched.
    @test d[2] ≈ 1.0
    @test d[3] ≈ 1.0
    # d_i > 1 exactly on tilted sectors.
    @test (d .> 1) == [true, false, false, true]
end

@testset "f1 tilt: cell_financing wires the tilt and rejects the stand-in" begin
    data = tiny_fixture()
    ψ = [0.7, 0.3]
    g = [0.07, 0.03]
    fin = cell_financing("F1", ψ, g, "tilt_g0_over_c0", data)
    @test fin isa PreferenceReallocation
    @test fin.shift ≈ f1_tilt_weights(data.household_baseline, ψ, g)
    err = try
        cell_financing("F1", ψ, g, "one_plus_psi", data)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("one_plus_psi", sprint(showerror, err))
    @test_throws ArgumentError cell_financing("F1", ψ, g, "bogus", data)
    @test_throws ArgumentError cell_financing("F9", ψ, g, "tilt_g0_over_c0", data)
end

@testset "f1 tilt: helper rejects an all-nonpositive baseline" begin
    @test_throws ArgumentError f1_tilt_weights([0.0, -1.0], [0.5, 0.5], [0.05, 0.05])
end

@testset "f1 tilt: helper validates inputs and floors the denominator" begin
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [0.7], [0.07, 0.03])
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [0.7, 0.3], [0.07])
    # Programme with zero mass over positive-baseline sectors.
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.0], [0.0, 1.0], [0.05, 0.05])
    # Negative ψ is rejected.
    @test_throws ArgumentError f1_tilt_weights([1.0, 1.0], [-0.1, 1.1], [0.1, 0.1])
    # Non-finite baseline is rejected (baseline may be zero/negative, but not NaN).
    @test_throws ArgumentError f1_tilt_weights([NaN, 1.0], [0.5, 0.5], [0.05, 0.05])
    @test_throws ArgumentError f1_tilt_weights([Inf, 1.0], [0.5, 0.5], [0.05, 0.05])
    # Non-finite or negative g is rejected.
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [0.5, 0.5], [NaN, 0.05])
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [0.5, 0.5], [Inf, 0.05])
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [0.5, 0.5], [-0.05, 0.05])
    # Non-finite ψ is rejected.
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [NaN, 0.5], [0.05, 0.05])
    @test_throws ArgumentError f1_tilt_weights([0.5, 0.5], [Inf, 0.5], [0.05, 0.05])
    # Denominator floor: tiny-but-positive baselines use max.(c0, 1e-12).
    d = f1_tilt_weights([1e-15, 1.0], [0.5, 0.5], [0.05, 0.05])
    @test d[1] ≈ 1 + 0.1 * 0.5 / 1e-12
    @test d[2] ≈ 1 + 0.1 * 0.5 / 1.0
end
