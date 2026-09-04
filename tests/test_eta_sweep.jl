using BeyondHulten
using Test

@testset "eta sweep result and accessors" begin
    data = tiny_fixture()
    shocks = Shocks([1.05, 1.0], ones(2), zeros(2))
    η_values = [0.0, 0.5]
    solutions = eta_sweep(data, shocks, 0.5, 0.5, 0.9, η_values; verbose=false)
    result = EtaSweepResult(Float64.(η_values), solutions)

    @test length(solutions) == length(η_values)
    @test result[1][1] == 0.0
    @test result[1][2] === solutions[1]
    @test length(real_gdp_sweep(result)) == 2
    @test length(nominal_gdp_sweep(result)) == 2
    @test size(sectoral_quantities(result)) == (2, 2)
    @test size(sectoral_prices(result)) == (2, 2)
    @test all(isfinite, real_gdp_sweep(result))

    reversed = eta_sweep(data, shocks, 0.5, 0.5, 0.9, reverse(η_values); verbose=false)
    by_eta = Dict(η => real_gdp_sweep(result)[i] for (i, η) in enumerate(η_values))
    @test all(isapprox(by_eta[η], real_gdp_sweep(EtaSweepResult(Float64.(reverse(η_values)), reversed))[i]; rtol=1e-8, atol=1e-10)
        for (i, η) in enumerate(reverse(η_values)))
end

@testset "eta sweep validation and contextual errors" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    @test_throws ArgumentError eta_sweep(data, shocks, 0.5, 0.5, 0.9, [Inf])
    @test_throws ArgumentError eta_sweep(data, shocks, 0.5, 0.5, 0.9, [51.])
    err = try
        eta_sweep(data, shocks, Inf, 0.5, 0.9, [0.])
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("θ", sprint(showerror, err))
end
