using BeyondHulten
using DataFrames
using LinearAlgebra
using Test

"""Deterministic 3-sector fixture with non-trivial conditional input shares.

Rows of `Ω` sum to 1 (user-by-supplier orientation). Accounting consistency:
`labor_share == λ .* factor_share` and `value_added == factor_share .* grossy`.
Uses the compact `Data` constructor; no disk reads.
"""
function three_sector_fixture()
    io = DataFrame("Sektoren" => ["a", "b", "c"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
    Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
    consumption_share = [0.5, 0.3, 0.2]
    factor_share = [0.6, 0.5, 0.4]
    λ = [1.2, 1.0, 0.8]
    labor_share = λ .* factor_share
    consumption_share_gross_output = [0.4, 0.35, 0.3]
    grossy = [10.0, 8.0, 6.0]
    value_added = factor_share .* grossy
    Data(io, Ω, consumption_share, factor_share, λ, labor_share,
        consumption_share_gross_output, grossy, value_added)
end

# Goldens re-pinned 2026-09-18 (ADR-0020 option C: sectoral wages at η = 0).
# Solved with default init (p=1, y=λ, w=1, F=0); no explicit `init` is passed,
# so these values also pin the default-init behavior. Tolerances are
# deliberately loose (atol=1e-5) so the tests only fail on genuine behavior
# change, not on solver-version noise.
#
# η = 0 is the SECTORAL-WAGE endpoint: the per-sector wage w_i is solved so that
# the frozen allocation is cost-minimizing at the equilibrium quantities. That
# replaces the retired single-wage pin, which could only REPORT the resulting
# factor-market gap. Two consequences are visible below: (i) the wage is a
# VECTOR, and it moves even in the unshocked cell, because this fixture's
# equilibrium quantities are not the baseline λ (q ≠ λ at both endpoints — a
# property of the fixture's compact calibration, not of the closure), so w_i
# must adjust to keep L^cm_i = L̄_i; (ii) the account closes exactly (canary gap
# ≈ 0 — see tests/test_external_closure.jl). η = 1 keeps its single wage and is
# untouched by ADR-0020.
const _GOLDEN_P_ETA0 = [1.065063296336466, 0.9556920527163134, 0.9133212699349667]
const _GOLDEN_Q_ETA0 = [1.2185742786684088, 0.9955914050894966, 0.7691557740796197]
const _GOLDEN_W_ETA0 = [1.0982897746921572, 0.9472841090297872, 0.844252211984721]
const _GOLDEN_F_ETA0 = 0.0
const _GOLDEN_RGDP_ETA0 = 0.996474972989093
const _GOLDEN_NGDP_ETA0 = 1.534571400128357

const _GOLDEN_P_ETA1 = [1.0, 1.0, 1.0]
const _GOLDEN_Q_ETA1 = [1.2822805578342904, 0.9651845775225595, 0.720098441345365]
const _GOLDEN_W_ETA1 = 1.0
const _GOLDEN_F_ETA1 = 0.0
const _GOLDEN_RGDP_ETA1 = 1.0
const _GOLDEN_NGDP_ETA1 = 1.54

# Legacy additive-shock compatibility path (autonomous + investment demand plus
# a sectoral supply shock), η = 0. Re-pinned 2026-09-18 (ADR-0020 option C):
# prices, quantities, the sectoral wages and real GDP all move, because the
# endpoint now solves the per-sector FOC instead of pinning one wage. F carries
# the legacy manna (the identity gap reads exactly p·(A+G) there and F absorbs
# it — see the manna testset in tests/test_external_closure.jl), so the canary
# gap is NOT zero in this cell by design.
const _GOLDEN_P_ADD = [0.9619860601261832, 1.0739502046103735, 0.9893540418786801]
const _GOLDEN_Q_ADD = [1.5075274008468775, 1.022770314723426, 0.791615546618296]
const _GOLDEN_W_ADD = [1.2651891823858499, 1.1234154023922183, 0.9687247326758176]
const _GOLDEN_F_ADD = -0.0988811763559975
const _GOLDEN_RGDP_ADD = 1.093347110093526

function _solve_mobile(data, shocks, η; kwargs...)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η;
        labor_bar=sum(data.labor_share), kwargs...)
    model, solve(model)
end

@testset "kernel regression: tiny fixture still solves" begin
    data = tiny_fixture()
    model, sol = _solve_mobile(data, Shocks(ones(2), ones(2), zeros(2)), 0.0)
    @test real_gdp(sol) ≈ 1 atol=1e-5
    @test nominal_gdp(sol) ≈ 1 atol=1e-5
    @test max_equilibrium_residual(sol) < 1e-5
end

@testset "kernel regression: mobile BF endpoints (η = 0 and 1)" begin
    data = three_sector_fixture()
    @test vec(sum(data.Ω_raw, dims=2)) ≈ ones(3) atol=1e-12
    @test data.labor_share ≈ data.λ .* data.factor_share atol=1e-12
    shocks = Shocks(ones(3), ones(3), zeros(3))

    model0, sol0 = _solve_mobile(data, shocks, 0.0)
    @test sol0.prices_raw ≈ _GOLDEN_P_ETA0 atol=1e-5
    @test sol0.quantities ≈ _GOLDEN_Q_ETA0 atol=1e-5
    @test sol0.wages_raw ≈ _GOLDEN_W_ETA0 atol=1e-5   # sectoral wages (ADR-0020)
    @test sol0.external_transfer ≈ _GOLDEN_F_ETA0 atol=1e-5
    @test real_gdp(sol0) ≈ _GOLDEN_RGDP_ETA0 atol=1e-5
    @test nominal_gdp(sol0) ≈ _GOLDEN_NGDP_ETA0 atol=1e-5

    model1, sol1 = _solve_mobile(data, shocks, 1.0)
    @test sol1.prices_raw ≈ _GOLDEN_P_ETA1 atol=1e-5
    @test sol1.quantities ≈ _GOLDEN_Q_ETA1 atol=1e-5
    @test sol1.wages_raw[1] ≈ _GOLDEN_W_ETA1 atol=1e-5
    @test sol1.external_transfer ≈ _GOLDEN_F_ETA1 atol=1e-5
    @test real_gdp(sol1) ≈ _GOLDEN_RGDP_ETA1 atol=1e-5
    @test nominal_gdp(sol1) ≈ _GOLDEN_NGDP_ETA1 atol=1e-5

    for (m, s) in ((model0, sol0), (model1, sol1))
        # η = 0 canonical vector is 3N+1 [p; y; w(1:N); F]; η = 1 is 2N+2.
        w = m.options.elasticities.η == 0.0 ? s.wages_raw : s.wages_raw[1]
        X = [s.prices_raw; s.quantities; w; s.external_transfer]
        @test maximum(abs, equilibrium_residuals(m, X)) < 1e-5
    end
    # η = 0 returns the frozen allocation by construction (the sectoral wages
    # are chosen so that it is cost-minimizing at the equilibrium quantities);
    # η = 1 reports the cost-minimizing demand at its single wage, which differs
    # from the frozen allocation — asserted below.
    @test sectoral_labor_demand(sol0.prices_raw, sol0.quantities, sol0.wages_raw, model0) ≈
        data.labor_share
    @test !isapprox(
        sectoral_labor_demand(sol1.prices_raw, sol1.quantities, sol1.wages_raw[1], model1),
        data.labor_share; atol=1e-6)
    @test_throws DomainError _solve_mobile(data, shocks, 0.5)
end

@testset "kernel regression: mobile with legacy additive shocks" begin
    data = three_sector_fixture()
    supply = ones(3)
    supply[1] = 1.2
    shocks = Shocks(supply, ones(3);
        autonomous_demand=[0.1, 0.0, 0.0],
        investment_shock=[0.0, 0.05, 0.0])
    model, sol = _solve_mobile(data, shocks, 0.0)
    @test sol.prices_raw ≈ _GOLDEN_P_ADD atol=1e-5
    @test sol.quantities ≈ _GOLDEN_Q_ADD atol=1e-5
    @test sol.wages_raw ≈ _GOLDEN_W_ADD atol=1e-5
    @test sol.external_transfer ≈ _GOLDEN_F_ADD atol=1e-5
    @test real_gdp(sol) ≈ _GOLDEN_RGDP_ADD atol=1e-5
    X = [sol.prices_raw; sol.quantities; sol.wages_raw; sol.external_transfer]
    @test maximum(abs, equilibrium_residuals(model, X)) < 1e-5
end

@testset "kernel regression: mobile contracts" begin
    data = three_sector_fixture()
    shocks = Shocks(ones(3), ones(3), zeros(3))
    model, sol = _solve_mobile(data, shocks, 1.0)
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
    @test maximum(abs, equilibrium_residuals(model, X)) < 1e-5

    labor = sectoral_labor_demand(sol.prices_raw, sol.quantities,
        sol.wages_raw[1], model)
    @test sum(sol.prices_raw .* sol.consumption) ≈
        sol.wages_raw[1] * sum(labor) atol=1e-9

    model_eta0 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0;
        labor_bar=sum(data.labor_share))
    @test sectoral_labor_demand(ones(3), ones(3), 1.0, model_eta0) ≈
        data.labor_share

    bf_model, bf_sol = _solve_mobile(data, shocks, 1.0)
    alpha_model, alpha_sol = _solve_mobile(data, shocks, 1.0)
    @test bf_sol.quantities ≈ alpha_sol.quantities atol=1e-10
    @test labor_closure(bf_model) isa FlexibleWageClosure
    @test labor_closure(model) isa FlexibleWageClosure
end

@testset "kernel regression: fixed GAMMA contracts only" begin
    # NOTE: exact fixed-closure numbers are intentionally NOT pinned. Phase 2
    # backports the frozen cbase2 formulation, which enforces all N
    # market-clearing equations and resolves the recorded "omitted N-th
    # market" open item — so fixed-closure levels are expected to move. Only
    # the closure contracts below must survive.
    data = three_sector_fixture()
    # Unanchored fixed-wage models must use η = 0: the η = 1 fixed system is
    # homogeneous and the scale-indeterminacy guard fires (tested below).
    model = mobile_labor_model(data, Shocks(ones(3), ones(3), zeros(3)),
        0.5, 0.5, 0.9, 0.0; closure=:fixed)
    sol = solve(model)
    @test sol.wages ≈ ones(3)
    @test maximum(abs,
        equilibrium_residuals(model, [sol.prices_raw; sol.quantities])) < 1e-5
    labor = sectoral_labor_demand(sol.prices_raw, sol.quantities, 1.0, model)
    @test sum(sol.prices_raw .* sol.consumption) ≈ sum(labor) atol=1e-9

    # ADR-0014/ADR-0017: manna does NOT anchor the η = 1 fixed system. Manna is
    # a constant added to final demand, so it leaves the actual clearing matrix
    # G untouched: on this closed fixture (m = s = 0) (I − G) is singular at
    # the baseline prices, 1'(I − G) = 0, and the solution set stays a line.
    # The retired heuristic admitted the case anyway and the solver returned a
    # point whose location was set by its path. The determinate counterpart
    # (same cell on a fixture with s > 0, where the actual clearing matrix
    # contracts) is asserted in tests/test_fixed_closure.jl.
    shocked = mobile_labor_model(data,
        Shocks(ones(3), ones(3); autonomous_demand=[0.1, 0.0, 0.0]),
        0.5, 0.5, 0.9, 1.0; closure=:fixed)
    err_manna = try
        solve(shocked)
        nothing
    catch e
        e
    end
    @test err_manna isa ArgumentError
    @test occursin("scale-indeterminate", sprint(showerror, err_manna))

    anchored = mobile_labor_model(data, Shocks(ones(3), ones(3), zeros(3)),
        0.5, 0.5, 0.9, 1.0; closure=:fixed)
    err = try
        solve(anchored)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("scale-indeterminate", sprint(showerror, err))
end
