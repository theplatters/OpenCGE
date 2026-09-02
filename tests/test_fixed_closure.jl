using BeyondHulten
using DataFrames
using Test

@testset "fixed wage closure" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    model = Model(data, shocks, MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, .5), 1., :fixed))
    @test labor_closure(model.options) isa FixedWageClosure
    @test labor_closure(model) isa FixedWageClosure
    sol = solve(model)
    @test real_gdp(sol) ≈ 1 atol=1e-5
    @test nominal_gdp(sol) ≈ 1 atol=1e-5
    @test sol.wages ≈ ones(2)
    @test maximum(abs, equilibrium_residuals(sol)) < 1e-5
    labor = sectoral_labor_demand(sol.prices_raw, sol.quantities, 1.0, model)
    @test sum(sol.prices_raw .* sol.consumption) ≈ sum(labor) atol=1e-10

    # Employment is an outcome, not a labor-market-clearing constraint, under
    # the fixed-wage closure.
    shocked = Model(data,
        Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        model.options)
    shocked_sol = solve(shocked)
    shocked_labor = sectoral_labor_demand(
        shocked_sol.prices_raw, shocked_sol.quantities, 1.0, shocked)
    @test maximum(abs, equilibrium_residuals(shocked_sol)) < 1e-5
    @test !isapprox(sum(shocked_labor), shocked.options.labor_bar; atol=1e-4)
    @test sum(shocked_sol.prices_raw .* shocked_sol.consumption) ≈
        sum(shocked_labor) atol=1e-9
end

@testset "Sobol result API" begin
    @test SobolResult === VarianceDecompositionResult
    grid = DataFrame(η=[0., 0., 1., 1.], ϵ=[0., 1., 0., 1.], θ=[0., 0., 0., 0.], σ=[0., 0., 0., 0.])
    result = BeyondHulten._compute_sobol_indices(["η", "ϵ", "θ", "σ"], grid,
        [0., 1., 1., 2.], "synthetic")
    @test result isa SobolResult
    @test result.S_f["η"] ≈ .5
    @test result.S_f["ϵ"] ≈ .5
    @test result.n_failed == 0

    public_result = variance_decomposition(
        tiny_fixture(),
        Shocks([1.1, 1.0], ones(2), zeros(2));
        η_values=[0.0, 1.0],
        ϵ_values=[0.5],
        θ_values=[0.5],
        σ_values=[0.9],
        verbose=false,
    )
    @test public_result isa SobolResult
    @test public_result.n_failed == 0
end
