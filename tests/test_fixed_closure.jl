using BeyondHulten
using DataFrames
using Test

@testset "fixed wage closure" begin
    data = tiny_fixture()
    shocks = Shocks(ones(2), ones(2), zeros(2))
    model = Model(data, shocks, MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 0.), 1., :fixed))
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
    # the fixed-wage closure. Use the mobile allocation (η = 1) so the reported
    # labor responds to the shock instead of being the constant η = 0 baseline
    # allocation.
    #
    # ADR-0014: the closed tiny_fixture has UNIT round-gain column sums
    # (A_bill/λ + (1−m)(1−s)fs = 1 exactly), so 1'(I − G) = 0 and its fixed
    # η = 1 system is scale-indeterminate — manna is a CONSTANT and cannot
    # remove a unit root, which is why the retired heuristic admitted this case
    # while the verified criterion rejects it. Give the fixture a positive
    # saving rate so the round-gain contracts (max colsum = 0.95 < 1) and the
    # system is determinate.
    d = data
    d_det = Data(d.io, d.Ω, d.Ω_raw, d.consumption_share, d.factor_share, d.λ,
        d.labor_share, d.consumption_share_gross_output, d.grossy, d.value_added,
        d.gross_output_basic, d.value_added_components, d.imports_intermediate,
        d.import_share, d.domestic_final_demand, d.gov_demand, d.household_baseline,
        d.import_margin, d.exo_demand, d.exports_demand, 0.1,
        d.A_bill, d.M_int, d.T_int, d.gdp_production, d.gdp_income, d.gdp_expenditure)
    shocked = Model(d_det,
        Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
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

@testset "fixed eta one: manna cannot anchor a closed fixture (ADR-0014)" begin
    # The retired heuristic admitted this case whenever manna was present. Manna
    # is a constant added to final demand: it leaves G untouched, so the unit
    # root survives (max round-gain column sum = 1.0 exactly) and the solution
    # set stays a line. The verified criterion rejects it.
    data = tiny_fixture()
    model = Model(data, Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    err = try
        solve(model)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("scale-indeterminate", sprint(showerror, err))
    # A positive saving rate contracts the round-gain and admits the same cell.
    d = data
    d_det = Data(d.io, d.Ω, d.Ω_raw, d.consumption_share, d.factor_share, d.λ,
        d.labor_share, d.consumption_share_gross_output, d.grossy, d.value_added,
        d.gross_output_basic, d.value_added_components, d.imports_intermediate,
        d.import_share, d.domestic_final_demand, d.gov_demand, d.household_baseline,
        d.import_margin, d.exo_demand, d.exports_demand, 0.1,
        d.A_bill, d.M_int, d.T_int, d.gdp_production, d.gdp_income, d.gdp_expenditure)
    ok = Model(d_det, Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    sol = solve(ok)
    @test maximum(abs, equilibrium_residuals(sol)) < 1e-5
end
