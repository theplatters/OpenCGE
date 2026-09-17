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

# Goldens captured on 2026-09-17 after the ADR-0010 port (N-1 clearing + CPI,
# residual external account) and the retirement of the allocation wedge with
# η ∈ {0, 1}. Solved with default init (p=1, y=λ, w=1); no explicit `init` is
# passed, so these values also pin the default-init behavior. Tolerances are
# deliberately loose (atol=1e-5) so the tests only fail on genuine behavior
# change, not on solver-version noise.
#
# η = 1 (fully mobile) solves the mobile system exactly. η = 0 (immobile
# baseline allocation) leaves the labor FOC unenforced and its solution carries
# the fixed-allocation gap in the omitted N-th market; the two endpoints
# therefore differ in general (they coincide only where the omitted market is
# redundant, e.g. the closed fixture at η = 1).
const _GOLDEN_P_ETA0 = [1.0, 1.0, 1.0]
const _GOLDEN_Q_ETA0 = [1.2965771781445423, 0.9799857844319968, 0.78266717964435]
const _GOLDEN_W_ETA0 = 1.0
const _GOLDEN_RGDP_ETA0 = 1.0
const _GOLDEN_NGDP_ETA0 = 1.54

const _GOLDEN_P_ETA1 = [1.0, 1.0, 1.0]
const _GOLDEN_Q_ETA1 = [1.2822805578342904, 0.9651845775225595, 0.720098441345365]
const _GOLDEN_W_ETA1 = 1.0
const _GOLDEN_RGDP_ETA1 = 1.0
const _GOLDEN_NGDP_ETA1 = 1.54

# Legacy additive-shock compatibility path (autonomous + investment demand plus
# a sectoral supply shock), η=0. Captured 2026-09-17.
const _GOLDEN_P_ADD = [0.90589586160246, 1.1101434216426687, 1.091897747460701]
const _GOLDEN_Q_ADD = [1.640761247830256, 1.0347653859537826, 0.6935523230378825]
const _GOLDEN_W_ADD = 1.159091891146103
const _GOLDEN_RGDP_ADD = 1.159091890922907

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
    @test sol0.wages_raw[1] ≈ _GOLDEN_W_ETA0 atol=1e-5
    @test real_gdp(sol0) ≈ _GOLDEN_RGDP_ETA0 atol=1e-5
    @test nominal_gdp(sol0) ≈ _GOLDEN_NGDP_ETA0 atol=1e-5

    model1, sol1 = _solve_mobile(data, shocks, 1.0)
    @test sol1.prices_raw ≈ _GOLDEN_P_ETA1 atol=1e-5
    @test sol1.quantities ≈ _GOLDEN_Q_ETA1 atol=1e-5
    @test sol1.wages_raw[1] ≈ _GOLDEN_W_ETA1 atol=1e-5
    @test real_gdp(sol1) ≈ _GOLDEN_RGDP_ETA1 atol=1e-5
    @test nominal_gdp(sol1) ≈ _GOLDEN_NGDP_ETA1 atol=1e-5

    for (m, s) in ((model0, sol0), (model1, sol1))
        X = [s.prices_raw; s.quantities; s.wages_raw[1]]
        @test maximum(abs, equilibrium_residuals(m, X)) < 1e-5
    end
    # η = 0 reports the baseline allocation; η = 1 reports the cost-minimizing
    # demand at the equilibrium (the two coincide at the baseline only). The
    # endpoints differ in aggregate because the omitted N-th market carries the
    # η = 0 fixed-allocation/factor-market gap (ADR-0010).
    @test sectoral_labor_demand(sol0.prices_raw, sol0.quantities, sol0.wages_raw[1], model0) ≈
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
    @test sol.wages_raw[1] ≈ _GOLDEN_W_ADD atol=1e-5
    @test real_gdp(sol) ≈ _GOLDEN_RGDP_ADD atol=1e-5
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
    @test maximum(abs, equilibrium_residuals(model, X)) < 1e-5
end

@testset "kernel regression: mobile contracts" begin
    data = three_sector_fixture()
    shocks = Shocks(ones(3), ones(3), zeros(3))
    model, sol = _solve_mobile(data, shocks, 1.0)
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
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

    # ADR-0014: manna does NOT anchor the η = 1 fixed system. Manna is a
    # constant added to final demand, so it leaves the round-gain matrix G
    # untouched: on this closed fixture (m = s = 0) the column sums are exactly
    # 1, 1'(I − G) = 0, and the solution set stays a line. The retired
    # heuristic admitted the case anyway and the solver returned a point whose
    # location was set by its path. The determinate counterpart (same cell on a
    # fixture with s > 0, where the round-gain contracts) is asserted in
    # tests/test_fixed_closure.jl.
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
