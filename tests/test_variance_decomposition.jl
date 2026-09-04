using BeyondHulten
using DataFrames
using Test

@testset "Sobol result API and balanced additive indices" begin
    @test SobolResult === VarianceDecompositionResult
    grid = DataFrame(η=[0., 0., 1., 1.], ϵ=[0., 1., 0., 1.])
    result = BeyondHulten._compute_sobol_indices(["η", "ϵ"], grid,
        [0., 2., 1., 3.], "synthetic")
    @test result isa SobolResult
    @test result.S_f["η"] ≈ .2
    @test result.S_f["ϵ"] ≈ .8
    @test result.ST_f["η"] ≈ .2
    @test result.ST_f["ϵ"] ≈ .8
    @test result.n_failed == 0
end

@testset "Sobol validation and incomplete grids" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    @test_throws ArgumentError variance_decomposition(data, shocks; output=:unknown, verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        η_values=Float64[], verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        η_values=[NaN], verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        η_values=[Inf], verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        η_values=[51.], verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        ϵ_values=Float64[], verbose=false)
    @test_throws ArgumentError variance_decomposition(data, shocks;
        θ_values=[NaN], verbose=false)
    @test_throws ArgumentError BeyondHulten._validate_grids(Float64[])
    @test_throws ArgumentError BeyondHulten._validate_grids([0., NaN])
    @test_throws ArgumentError BeyondHulten._validate_grids([Inf]; eta=true)
    @test_throws ArgumentError BeyondHulten._validate_grids([51.]; eta=true)

    grid = DataFrame(η=[0., 0., 1., 1.], ϵ=[0., 1., 0., 1.])
    incomplete = BeyondHulten._compute_sobol_indices(["η", "ϵ"], grid,
        [0., 2., 1., NaN], "incomplete")
    @test incomplete.n_failed == 1
    # The expected total-order values use the actual three valid observations:
    # grand mean is 1, SS_total = 2; SS_except_η = 2*(1/2)^2 + 1 = 3/2.
    @test incomplete.ST_f["η"] ≈ 1 - (3 / 2) / 2
    # At fixed η the valid group means are 1 and 1, so SS_except_ϵ = 0.
    @test incomplete.ST_f["ϵ"] ≈ 1

    invalid = BeyondHulten._compute_sobol_indices(["η", "ϵ"], grid,
        [NaN, Inf, NaN, Inf], "invalid")
    @test invalid.n_failed == 4
    @test_nowarn summary_table(invalid)
end

@testset "public variance decomposition" begin
    result = variance_decomposition(tiny_fixture(),
        Shocks([1.1, 1.0], ones(2), zeros(2));
        η_values=[0.0, 1.0], ϵ_values=[0.5], θ_values=[0.5],
        σ_values=[0.9], verbose=false)
    @test result isa SobolResult
    @test result.n_failed == 0
    @test result.output_name == "real_gdp"
end
