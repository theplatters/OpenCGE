using BeyondHulten
using DataFrames
using Test

@testset "mobile labor closure and reallocation" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    model = mobile_labor_model(data, shocks, .5, .5, .9, 0.; labor_bar=1.)
    @test labor_closure(model.options) isa FlexibleWageClosure
    @test labor_closure(model) isa FlexibleWageClosure
    sol = solve(model)
    @test real_gdp(sol) ≈ 1 atol=1e-5
    @test nominal_gdp(sol) ≈ 1 atol=1e-5
    @test max_equilibrium_residual(sol) < 1e-5
    labor = sectoral_labor_demand(sol.prices_raw, sol.quantities,
        sol.wages_raw[1], model)
    @test sum(sol.prices_raw .* sol.consumption) ≈
        sol.wages_raw[1] * sum(labor) atol=1e-10

    e = model.options.elasticities
    @test sectoral_labor_demand(ones(2), ones(2), 1., Model(data, shocks,
        MobileLaborCES(e, 1., :mobile))) ≈ data.labor_share
    m0 = Model(data, shocks, MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 0.), 1., :mobile))
    m1 = Model(data, shocks, MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :mobile))
    p, y = [1.2, .8], [1.1, .9]
    @test sectoral_labor_demand(p, y, 1., m0) ≈ data.labor_share
    @test sectoral_labor_demand(p, y, 1., m1) ≈
        (data.factor_share .* y .* sqrt.(p)) atol=1e-12
    @test MobileLaborCES(e, 1., FlexibleWageClosure()).closure == :mobile
    @test mobile_labor_model(data, shocks, .5, .5, .9, 0.;
        closure=FixedWageClosure()).options.closure == :fixed
    @test mobile_labor_model(data, shocks, .5, .5, .9, 0.;
        closure=:fixed).options.labor_bar == sum(data.labor_share)
    @test_throws ArgumentError mobile_labor_model(data, shocks, .5, .5, .9, 0.;
        labor_bar=1., closure=:fixed)
    @test_throws ArgumentError mobile_labor_model(data, shocks, .5, .5, .9, 0.;
        labor_bar=1., closure=FixedWageClosure())
    @test_throws ArgumentError MobileLaborCES(e, 1., :unknown)
    @test_throws ArgumentError MobileLaborCES(e, 1.; closure="mobile")
    failure = try
        solve(model; init=fill(NaN, 5))
        nothing
    catch err
        err
    end
    @test failure isa ErrorException
    @test occursin("did not converge", sprint(showerror, failure))
end

@testset "mobile labor numerical safeguards and eta endpoints" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    model = Model(data, shocks,
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 0.), 1., :mobile))

    # η = 0 returns the baseline allocation; η = 1 the cost-minimizing demand
    # (equal to the baseline at p = y = w = 1). Intermediate values are
    # retired (ADR-0010) and rejected loudly.
    @test sectoral_labor_demand(ones(2), ones(2), 1., model) ≈ data.labor_share
    m1 = Model(data, shocks,
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :mobile))
    @test sectoral_labor_demand(ones(2), ones(2), 1., m1) ≈ data.labor_share
    @test_throws ArgumentError sectoral_labor_demand(ones(2), ones(2), 1.,
        Model(data, shocks, MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, Inf), 1., :mobile)))
    for retired in (0.5, 50.)
        @test_throws DomainError sectoral_labor_demand(ones(2), ones(2), 1.,
            Model(data, shocks,
                MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, retired), 1., :mobile)))
    end
    @test economy_wide_wage(ones(2), ones(2), data.labor_share, model) ≈ 1.
end

@testset "legacy closure taxonomy and compatibility" begin
    callback = m -> m.data.labor_share
    options = CES(CESElasticities(.5, .5, .9), callback, false)
    @test options.labor_slack === callback
    @test labor_closure(options) isa ExogenousLaborClosure
    legacy_model = Model(tiny_fixture(), Shocks(ones(2), ones(2), zeros(2)), options)
    @test labor_closure(legacy_model) isa ExogenousLaborClosure
    @test labor_closure(options).callback === callback
    cd = CobbDouglas(CobbDouglasElasticities([0.5], [0.5]), callback)
    @test labor_closure(cd).callback === callback
    @test CES(CESElasticities(.5, .5, .9); labor_slack=callback).labor_slack === callback
    positional = Shocks(ones(2), ones(2), zeros(2))
    @test positional.autonomous_demand == zeros(2)
    @test Shocks(ones(2), ones(2), zeros(2), [.1, .2], zeros(2)).autonomous_demand == [.1, .2]
    @test Shocks(ones(2), ones(2); autonomous_demand=[.1, .2]).autonomous_demand == [.1, .2]
    @test_throws DimensionMismatch Shocks(ones(2), ones(3), zeros(2))
end
