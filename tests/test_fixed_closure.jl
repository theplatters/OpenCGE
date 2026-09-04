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
    @test max_equilibrium_residual(sol) < 1e-5
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
    @test nominal_gdp(shocked_sol) ≈ sum(shocked_labor) atol=1e-9
    @test sum(shocked_sol.prices_raw .* shocked_sol.consumption) ≈
        sum(shocked_labor) atol=1e-9
end

@testset "fixed eta near one scale validation" begin
    data = tiny_fixture()
    model = Model(data, Shocks(ones(2), ones(2), zeros(2)),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, .999999), 1., :fixed))
    err = try
        solve(model)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("scale-indeterminate", sprint(showerror, err))
end

@testset "fixed eta one anchor validation" begin
    data = tiny_fixture()
    model = Model(data, Shocks(ones(2), ones(2), zeros(2)),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    err = try
        solve(model)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("autonomous or investment", sprint(showerror, err))
end
