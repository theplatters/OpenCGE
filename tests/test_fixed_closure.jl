using BeyondHulten
using DataFrames
using LinearAlgebra
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
    # ADR-0014/ADR-0017: the closed tiny_fixture has a singular (I − G) at the
    # baseline prices (the actual clearing matrix G has column sums exactly 1,
    # so 1'(I − G) = 0), hence its fixed η = 1 system is scale-indeterminate —
    # manna is a CONSTANT and cannot remove a unit root, which is why the
    # retired heuristic admitted this case while the verified criterion
    # rejects it. Give the fixture a positive saving rate so the actual
    # clearing matrix contracts (max colsum = 0.95 < 1) and the system is
    # determinate.
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
    # Household expenditure is (1 − s) × income; the fixture carries s = 0.1 so
    # the actual clearing matrix contracts (ADR-0014/ADR-0017). With s = 0 the
    # factor would be 1.
    @test sum(shocked_sol.prices_raw .* shocked_sol.consumption) ≈
        (1 - 0.1) * sum(shocked_labor) atol=1e-9
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
    # ADR-0017: the closed tiny_fixture is rejected because the ACTUAL clearing
    # matrix (I − G) is singular at the baseline prices (a unit root, hence a
    # continuum of solutions). The guard diagnoses the matrix, not the absent
    # anchor: additive demand is a constant and cannot remove the singularity.
    data = tiny_fixture()
    model = Model(data, Shocks(ones(2), ones(2), zeros(2)),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    err = try
        solve(model)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("scale-indeterminate", sprint(showerror, err))
    @test occursin("singular", sprint(showerror, err))
end

@testset "fixed eta one: manna cannot anchor a closed fixture (ADR-0014)" begin
    # The retired heuristic admitted this case whenever manna was present. Manna
    # is a constant added to final demand: it leaves G untouched, so the unit
    # root survives (the actual (I − G) is singular at the baseline prices)
    # and the solution set stays a line. The verified criterion rejects it.
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
    # A positive saving rate contracts the actual clearing matrix and admits
    # the same cell.
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

@testset "fixed eta one: unit column sum without a unit root (ADR-0017)" begin
    # A maximum column sum of 1 does NOT imply a unit root. With
    # A_bill = [0.5, 0.4], M_int = [0.0, 0.1] (s = m = 0) the actual clearing
    # matrix has column sums [1.0, 0.9] but det(I − G) = 0.025 ≠ 0, so the
    # fixed η = 1 system is determinate and solves exactly at p = [1, 1],
    # y = [0.7, 0.5]. The retired closed-form trigger rejected this cell.
    d = tiny_fixture()
    vals = Any[getfield(d, f) for f in fieldnames(Data)]
    vals[findfirst(==(:A_bill), fieldnames(Data))] = [0.5, 0.4]
    vals[findfirst(==(:M_int), fieldnames(Data))] = [0.0, 0.1]
    d2 = Data(vals...)
    model = Model(d2, Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    G, _ = BeyondHulten._fixed_clearing_affine(model, ones(2))
    @test vec(sum(G; dims=1)) ≈ [1.0, 0.9] atol=1e-12
    @test det(I - G) ≈ 0.025 atol=1e-12
    sol = solve(model)
    @test maximum(abs,
        equilibrium_residuals(model, [sol.prices_raw; sol.quantities])) < 1e-6
    @test sol.prices_raw ≈ ones(2) atol=1e-9
    @test sol.quantities ≈ [0.7, 0.5] atol=1e-9
end

@testset "fixed eta one: heterogeneous import margins enter at the spending sector (ADR-0017)" begin
    # Import margins enter at the SPENDING sector through the household
    # expenditure composition, so the actual column sums are [0.95, 0.95] —
    # not the retired closed form's [1.0, 0.9], which applied each sector's
    # own margin. The retired trigger read max = 1.0 and rejected this cell;
    # the actual matrix contracts and the cell solves at y = [0.6, 0.4].
    d = tiny_fixture()
    vals = Any[getfield(d, f) for f in fieldnames(Data)]
    vals[findfirst(==(:import_margin), fieldnames(Data))] = [0.0, 0.2]
    d2 = Data(vals...)
    model = Model(d2, Shocks(ones(2), ones(2); autonomous_demand=[0.1, 0.0]),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    G, _ = BeyondHulten._fixed_clearing_affine(model, ones(2))
    @test vec(sum(G; dims=1)) ≈ [0.95, 0.95] atol=1e-12
    old = d2.A_bill ./ d2.λ .+ (1.0 .- d2.import_margin) .*
        (1.0 - d2.saving_rate) .* d2.factor_share
    @test old ≈ [1.0, 0.9] atol=1e-12
    sol = solve(model)
    @test maximum(abs,
        equilibrium_residuals(model, [sol.prices_raw; sol.quantities])) < 1e-6
    @test sol.quantities ≈ [0.6, 0.4] atol=1e-9
end

@testset "fixed eta one: a determinate but non-positive candidate is rejected (ADR-0017)" begin
    # The testset-1 calibration without manna (no exogenous demand at all) has
    # a NONSINGULAR (I − G) — det = 0.025 — but the unique clearing candidate
    # is the degenerate y = 0, which no equilibrium can realise. The guard
    # takes the positivity branch (not the singularity branch).
    d = tiny_fixture()
    vals = Any[getfield(d, f) for f in fieldnames(Data)]
    vals[findfirst(==(:A_bill), fieldnames(Data))] = [0.5, 0.4]
    vals[findfirst(==(:M_int), fieldnames(Data))] = [0.0, 0.1]
    d2 = Data(vals...)
    model = Model(d2, Shocks(ones(2), ones(2), zeros(2)),
        MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :fixed))
    G, c = BeyondHulten._fixed_clearing_affine(model, ones(2))
    @test abs(det(I - G)) > 1e-6
    @test maximum(abs, (I - G) \ c) < 1e-9
    err = try
        solve(model)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("positive", sprint(showerror, err))
end
