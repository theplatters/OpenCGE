using BeyondHulten
using DataFrames
using LinearAlgebra
using Test

# Contract tests for the Phase 2 promoted closures (ADR-0005): BETA/DELTA and
# the F1/F2/F3 financing closures. These are the evidence for upgrading
# BETA/DELTA/F1-F3 from `implemented` to `tested` in registry/closures.toml.
# Model: the assertions in cbase2/scripts/verify_v3.jl (frozen, read-only).
# All solves are deterministic 3-sector systems with default init unless a
# warm start is noted; no `read_data`, no disk writes.

"""Small v3-style synthetic fixture with an internally consistent open-absorption calibration.

Recipe: baseline wage income normalised to 1 (`sum(labor_share) == 1`);
uniform import margin `m`; small government/investment/export blocks;
`c0_dom = λ − Ω'((1−fs)λ) − (1−m)(gG+inv) − expo` (clamped positive by
construction values); `s`, `ω`, `household_baseline`, `consumption_share`
derived from it. Satisfies the finiteness gate columnwise.
"""
function v3_fixture()
    io = DataFrame("Sektoren" => ["a", "b", "c"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
    Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
    factor_share = [0.6, 0.5, 0.4]
    grossy = [2.0, 1.5, 1.0]
    value_added = factor_share .* grossy
    gdp = sum(value_added)
    λ = grossy ./ gdp
    labor_share = λ .* factor_share          # sums to 1 by construction
    m = 0.2
    gov = [0.03, 0.015, 0.005]               # Σ ≈ 0.05
    inv = [0.03, 0.075, 0.045]               # Σ = 0.15
    expo = [0.075, 0.045, 0.03]              # Σ = 0.15
    c0_dom = λ .- Ω' * ((1 .- factor_share) .* λ) .-
        (1 - m) .* (gov .+ inv) .- expo
    c0_gross = c0_dom ./ (1 - m)
    saving_rate = 1 - sum(c0_gross) / (1 - sum(gov))
    ω = c0_gross ./ sum(c0_gross)
    data = Data(io, Ω, Ω, ω, factor_share, λ, labor_share, ω, grossy,
        value_added, grossy, DataFrame(), zeros(3), zeros(3), zeros(3),
        gov, c0_gross, fill(m, 3), inv, expo, saving_rate,
        (1 .- factor_share) .* λ, zeros(3), zeros(3), gdp, gdp, gdp)
    # Small positive programme bundle on sector 1; strictly positive F1 shift.
    (; data = data, g = [0.02, 0.0, 0.0], shift = [1.5, 1.0, 0.8],
        saving_rate = saving_rate)
end

const _V3_θ, _V3_ϵ, _V3_σ, _V3_η = 1.0, 0.5, 0.9, 1.0
_v3_shocks() = Shocks(ones(3), ones(3), zeros(3))

@testset "promoted closures: v3-style fixture consistency" begin
    fx = v3_fixture()
    (; data) = fx
    @test sum(data.labor_share) ≈ 1.0 atol=1e-12
    @test data.labor_share ≈ data.λ .* data.factor_share atol=1e-12
    @test all(>(0), data.household_baseline)
    @test data.consumption_share ≈
        data.household_baseline ./ sum(data.household_baseline) atol=1e-12
    @test 0 < data.saving_rate < 1
    # Finiteness gate (task recipe): columnwise below 1 ...
    @test all((1 .- data.factor_share) .+
        (1 .- data.import_margin) .* (1 - data.saving_rate) .*
        data.factor_share .< 1)
    # ... and the true gain-matrix gate inside `leontief_multiplier` holds
    # (it errors when column sums reach 1).
    for mode in (:F2, :F3)
        ana = leontief_multiplier(data, fx.g; mode = mode)
        @test all(isfinite, ana.y) && all(isfinite, [ana.L, ana.F])
    end
end

@testset "promoted closures: BETA validation and ALPHA nesting" begin
    # `ElasticLaborClosure` validation is unchanged from cbase2: η_s finite
    # and non-negative, w0 positive.
    @test_throws ArgumentError ElasticLaborClosure(-0.5)
    @test_throws ArgumentError ElasticLaborClosure(Inf)
    @test_throws ArgumentError ElasticLaborClosure(NaN)
    @test_throws ArgumentError ElasticLaborClosure(0.5, 0.0)
    @test_throws ArgumentError ElasticLaborClosure(0.5, -1.0)
    @test ElasticLaborClosure(0.5) == ElasticLaborClosure(0.5, 1.0)

    fx = v3_fixture()
    shocks = _v3_shocks()
    beta = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        eta_s = 0.5)
    @test labor_closure(beta) isa ElasticLaborClosure
    @test labor_closure(beta.options) isa ElasticLaborClosure
    @test labor_closure(beta.options).η_s ≈ 0.5
    @test labor_closure(beta_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
        _V3_η; eta_s = 0.5)) isa ElasticLaborClosure
    # An explicit supply elasticity forces the :beta closure ...
    @test beta.options.closure === :beta
    # ... and `:beta` is accepted as a closure symbol directly.
    @test labor_closure(mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ,
        _V3_σ, _V3_η; closure = :beta)) isa ElasticLaborClosure

    # η_s = 0 reproduces the ALPHA/:mobile solution exactly (registry
    # formulation: vertical supply is the η_s = 0 limit of the BETA curve).
    fin = ExternalDebt(fx.g)
    m0 = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin, eta_s = 0.0)
    s0 = solve(m0)
    ma = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    sa = solve(ma)
    @test s0.prices_raw ≈ sa.prices_raw atol=1e-8
    @test s0.quantities ≈ sa.quantities atol=1e-8
    @test s0.wages_raw[1] ≈ sa.wages_raw[1] atol=1e-8
end

@testset "promoted closures: BETA supply curve (verify_v3 §4)" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    fin = ExternalDebt(fx.g)
    # Warm start for the continuation, as in verify_v3 (deterministic: the
    # linear fixed point is not needed; the ALPHA-regime solve is enough).
    ref = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    ref_sol = solve(ref)
    init = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]; ref_sol.external_transfer]
    lbar = sum(fx.data.labor_share)
    for η_s in (0.5, 1.0)
        sol = solve_beta(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
            financing = fin, eta_s = η_s, init = init)
        mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
            _V3_η; financing = fin, eta_s = η_s)
        w = sol.wages_raw[1]
        X = [sol.prices_raw; sol.quantities; w; sol.external_transfer]
        @test maximum(abs, equilibrium_residuals(mdl, X)) < 1e-6
        L = sum(sectoral_labor_demand(sol.prices_raw, sol.quantities, w, mdl))
        # Registry formulation Σ L_i = L̄·((w/P)/(w0/P0))^η_s with the CPI
        # numeraire (P = 1) and the w0 = 1 anchor. This is the enforced
        # labour-market equation; the former log(L/lbar)/log(w) "implied
        # elasticity" was just this residual solved for η_s (circular, review
        # §3.4) and is degenerate whenever the equilibrium keeps w = 1.
        @test L ≈ lbar * (w / 1.0)^η_s atol=1e-8
    end
end

@testset "promoted closures: DELTA equivalence (verify_v3 §3)" begin
    # `delta_elasticities(ε)` is the Leontief limit (θ, ϵ, σ → 0⁺) at η = 1.
    @test (e = delta_elasticities(1e-4);
        e.θ ≈ 1e-4 && e.ϵ ≈ 1e-4 && e.σ ≈ 1e-4 && e.η ≈ 1.0)
    fx = v3_fixture()
    shocks = _v3_shocks()   # demand-only shocks (A = 1)
    for mode in (:F2, :F3)
        fin = mode === :F2 ? TaxFinanced(fx.g) : ExternalDebt(fx.g)
        ana = leontief_multiplier(fx.data, fx.g; mode = mode)
        mdl = delta_model(fx.data, shocks; ε = 1e-4, financing = fin)
        sol = solve(mdl)
        rmax = maximum(abs, equilibrium_residuals(mdl,
            [sol.prices_raw; sol.quantities]))
        @test rmax < 1e-6
        rel = maximum(abs.(sol.quantities .- ana.y)) / maximum(abs.(ana.y))
        @test rel < 5e-3
        # ε-convergence toward the analytic system as ε → 0.
        sol3 = solve(delta_model(fx.data, shocks; ε = 1e-3, financing = fin))
        e3 = maximum(abs.(sol3.quantities .- ana.y))
        e4 = maximum(abs.(sol.quantities .- ana.y))
        @test e4 <= e3
        # Fixed-price multiplier structure: demand-only shocks pin p = 1.
        @test maximum(abs.(sol.prices_raw .- 1)) < 1e-6
    end
end

@testset "promoted closures: F1 compositional budget" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    fin = PreferenceReallocation(fx.shift)
    # Registry formulation β̃ᵢ = βᵢdᵢ/Σβⱼdⱼ side: the hook multiplies the
    # shift into the demand weights, others pass through unchanged.
    @test preference_weights(fin, shocks.demand_shock) ≈
        shocks.demand_shock .* fx.shift
    @test preference_weights(NoFinancing(), shocks.demand_shock) ≈
        shocks.demand_shock
    mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    sol = solve(mdl)
    p, w = sol.prices_raw, sol.wages_raw[1]
    @test maximum(abs, equilibrium_residuals(mdl,
        [p; sol.quantities; w; sol.external_transfer])) < 1e-5
    L = sum(sectoral_labor_demand(p, sol.quantities, w, mdl))
    E = household_expenditure(fin, mdl, w * L, p, L;
        external_transfer = sol.external_transfer)
    # Registry formulation Σ pᵢcᵢʰ = Eʰ (gross household consumption is
    # (1−s)E with the v3 saving leak; verify_v3 headline identity).
    @test dot(p, sol.consumption) ≈ (1 - fx.saving_rate) * E atol=1e-9
    # Compositional shift: expenditure shares tilt toward high-shift sectors
    # relative to the composition-neutral stand-in.
    neutral = solve(mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
        _V3_η; financing = PreferenceReallocation(ones(3))))
    sh = p .* sol.consumption ./ sum(p .* sol.consumption)
    sh0 = neutral.prices_raw .* neutral.consumption ./
        sum(neutral.prices_raw .* neutral.consumption)
    @test sh[1] > sh0[1]
    @test sh[3] < sh0[3]
    # F1 carries no additive demand (it must not anchor fixed-wage scale).
    @test additive_demand(fin, 3) ≈ zeros(3)
end

@testset "promoted closures: F2 tax-financed budget" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    fin = TaxFinanced(fx.g)
    mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    sol = solve(mdl)
    p, w = sol.prices_raw, sol.wages_raw[1]
    @test maximum(abs, equilibrium_residuals(mdl,
        [p; sol.quantities; w; sol.external_transfer])) < 1e-5
    L = sum(sectoral_labor_demand(p, sol.quantities, w, mdl))
    E = household_expenditure(fin, mdl, w * L, p, L;
        external_transfer = sol.external_transfer)
    # Registry formulation Σ pᵢgᵢ = T(p): the balanced-budget rule makes the
    # household pay baseline government plus the programme at current prices.
    # The external transfer enters AFTER tax (ADR-0019, decision D1), so the
    # expenditure identity carries + F.
    @test E ≈ w * L - (dot(p, fx.data.gov_demand) + dot(p, fx.g)) +
        sol.external_transfer atol=1e-9
    @test dot(p, sol.consumption) ≈ (1 - fx.saving_rate) * E atol=1e-9
    @test tau_rate(fin, mdl, p, w, L) > 0
    # `public_budget` reports baseline gG plus the programme at current
    # prices (its implemented definition; review §2.7 caveat is recorded in
    # the registry notes, not re-litigated here).
    @test public_budget(fin, mdl, p) ≈
        sum(fx.data.gov_demand) + dot(fx.g, p) atol=1e-9
    @test additive_demand(fin, 3) ≈ fx.g
end

@testset "promoted closures: F3 external budget" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    fin = ExternalDebt(fx.g)
    mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    sol = solve(mdl)
    p, w = sol.prices_raw, sol.wages_raw[1]
    @test maximum(abs, equilibrium_residuals(mdl,
        [p; sol.quantities; w; sol.external_transfer])) < 1e-5
    L = sum(sectoral_labor_demand(p, sol.quantities, w, mdl))
    E = household_expenditure(fin, mdl, w * L, p, L;
        external_transfer = sol.external_transfer)
    # Registry formulation: the programme is externally financed, so the
    # household is untaxed for it (only baseline gG is levied). The mobile
    # equilibrium transfer (ADR-0019) enters expenditure after tax, hence + F.
    @test E ≈ w * L - dot(p, fx.data.gov_demand) + sol.external_transfer atol=1e-9
    @test dot(p, sol.consumption) ≈ (1 - fx.saving_rate) * E atol=1e-9
    # `external_balance` records the programme's import content Σ p·m·g
    # (its implemented definition; review §2.8 caveat is recorded in the
    # registry notes, not re-litigated here).
    @test external_balance(fin, mdl, p) ≈
        dot(p, fx.data.import_margin .* fx.g) atol=1e-12
    @test additive_demand(fin, 3) ≈ fx.g
    @test additive_demand(NoFinancing(), 3) ≈ zeros(3)
end

@testset "promoted closures: external-account canary (review 2.1)" begin
    # ADR-0019: every regime enforces ALL N goods-market clearings with the
    # explicit external transfer F threaded into demand. At a mobile (η = 1)
    # equilibrium every clearing residual is ~0, so `dot(p, cl) ≈ 0` and the
    # canary identity gap `diff = S + T + M − (I+X) − (F + B_gov)` is ≈ 0. At
    # η = 0 the clearings still hold (with the F = 0 pin), but the canary
    # retains the documented fixed-allocation/factor-market gap (zero-profit
    # prices the cost-minimizing labour demand, not the frozen baseline
    # allocation): reported, not gated.
    fx = v3_fixture()
    shocks = _v3_shocks()
    for fin in (NoFinancing(), TaxFinanced(fx.g), ExternalDebt(fx.g))
        for η in (0.0, 1.0)
            mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, η;
                financing = fin)
            sol = solve(mdl)
            X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
            @test maximum(abs, equilibrium_residuals(mdl, X)) < 1e-6
            cl = market_clearing_residuals(mdl, X)
            can = external_balance_canary(mdl, X)
            # All N clearings are enforced at both endpoints.
            @test maximum(abs, cl) < 1e-6
            if η == 1.0
                @test dot(sol.prices_raw, cl) ≈ 0 atol=1e-9
                @test can.diff ≈ 0 atol=1e-9
            else
                # The η = 0 factor-market gap is a reported quantity, not a
                # gate: it must be finite (its magnitude is pinned on the real
                # table in tests/test_external_closure.jl).
                @test isfinite(can.diff)
            end
        end
    end

    # ADR-0013: the intermediate-tax leak (row 75) enters the canary linearly.
    # T_int is not part of the demand system, so the solution is unchanged; the
    # canary's diff must move by exactly the valued tax term.
    #
    # A synthetic fixture cannot demonstrate the *closure* of the identity: its
    # calibration is consistent under ANY split of the purchaser-price bill
    # (only A + M + T = (1−fs)·λ matters), whereas on the real table rows 73/74/
    # 75 fix the split and leaving row 75 out breaks the identity by exactly it.
    # The end-to-end contract is therefore asserted on the real table in
    # tests/test_calibration.jl ("ADR-0013 contract").
    d = fx.data
    T = fill(0.01, 3)
    d_t = Data(d.io, d.Ω, d.Ω_raw, d.consumption_share, d.factor_share, d.λ,
        d.labor_share, d.consumption_share_gross_output, d.grossy, d.value_added,
        d.gross_output_basic, d.value_added_components, d.imports_intermediate,
        d.import_share, d.domestic_final_demand, d.gov_demand, d.household_baseline,
        d.import_margin, d.exo_demand, d.exports_demand, d.saving_rate,
        d.A_bill, d.M_int, T, d.gdp_production, d.gdp_income, d.gdp_expenditure)
    mdl = mobile_labor_model(d, shocks, _V3_θ, _V3_ϵ, _V3_σ, 1.0)
    mdl_t = mobile_labor_model(d_t, shocks, _V3_θ, _V3_ϵ, _V3_σ, 1.0)
    sol = solve(mdl)
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
    can = external_balance_canary(mdl, X)
    can_t = external_balance_canary(mdl_t, X)
    @test can.T ≈ 0.0 atol = 1e-12
    @test can_t.T ≈ dot(sol.prices_raw .* (T ./ d.λ), sol.quantities) atol = 1e-12
    @test can_t.diff - can.diff ≈ can_t.T atol = 1e-12
    # The demand system is untouched by T_int, so the solution — quantities
    # and the external transfer — is identical.
    sol_t = solve(mdl_t)
    @test sol_t.quantities ≈ sol.quantities atol = 1e-12
    @test sol_t.external_transfer ≈ sol.external_transfer atol = 1e-12
end

@testset "promoted closures: intermediate-leak valuation is CES-consistent (ADR-0016)" begin
    # One-sector exact root: fs = 0.5, λ = 2, A_bill + M_int + T_int = 1.0 =
    # (1−fs)·λ, supply shock a = 4, ϵ = θ = 0.5 so the CES bill factor is
    # k = p^ϵ · a^(ϵ−1) · P^(1-ϵ) = 0.5. The imported intermediates and the
    # intermediate taxes are shares of the SAME bundle the domestic bill is
    # charged from, so both leaks must scale with k. Before ADR-0016 they were
    # valued at p alone and the canary read 3.0 instead of the exact 1.5.
    io = DataFrame("Sektoren" => ["a"], "Letzte Verwendung von Gütern zusammen" => [1.0])
    d0 = Data(io, ones(1, 1), [1.0], [0.5], [2.0], [1.0], [1.0], [2.0], [1.0])
    vals = Any[getfield(d0, f) for f in fieldnames(Data)]
    vals[findfirst(==(:A_bill), fieldnames(Data))] = [0.5]
    vals[findfirst(==(:M_int), fieldnames(Data))] = [0.25]
    vals[findfirst(==(:T_int), fieldnames(Data))] = [0.25]
    d = Data(vals...)
    mdl = mobile_labor_model(d, Shocks([4.0], [1.0], [0.0]), 0.5, 0.5, 0.9, 1.0)
    # Legacy 2N+1 vector (F = 0): zero-profit, labour and CPI hold exactly;
    # the single clearing residual is the unfinanced 1.5, and the canary reads
    # the same gap. Booking F = 1.5 closes the account exactly (ADR-0019).
    X = [1.0, 12.0, 9.0]
    r = equilibrium_residuals(mdl, X)
    @test r[1] ≈ 0 atol = 1e-12
    @test r[2] ≈ 1.5 atol = 1e-12
    @test r[3] ≈ 0 atol = 1e-12
    @test r[4] ≈ 0 atol = 1e-12
    @test dot([1.0], market_clearing_residuals(mdl, X)) ≈ 1.5 atol = 1e-12
    @test external_balance_canary(mdl, X).diff ≈ 1.5 atol = 1e-12
    @test maximum(abs, equilibrium_residuals(mdl, [1.0, 12.0, 9.0, 1.5])) < 1e-12
end

@testset "promoted closures: fixed-wage financing anchor at η = 1" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    N = length(fx.data.factor_share)
    # F2/F3 carry an additive anchor, so the fixed η = 1 system solves ...
    for fin in (TaxFinanced(fx.g), ExternalDebt(fx.g))
        mdl = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, 1.0;
            closure = :fixed, financing = fin)
        sol = solve(mdl)
        r = equilibrium_residuals(mdl, [sol.prices_raw; sol.quantities])
        # ... with all N market clearings enforced (length 2N residual).
        @test length(r) == 2N
        @test maximum(abs, r) < 1e-6
    end
    # ... while unanchored demand still throws the legacy guard.
    bare = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, 1.0;
        closure = :fixed)
    # ADR-0014/ADR-0017: the v3 fixture is an OPEN calibration (m = 0.2,
    # s ≈ 0.092), so the actual clearing matrix contracts (max column sum
    # 0.8905 < 1) and the fixed η = 1 system is DETERMINATE — it solves with
    # no additive anchor. The retired heuristic rejected it (while admitting
    # CLOSED fixtures whenever manna was present); the verified criterion does
    # the opposite, and the rejections are pinned in
    # tests/test_fixed_closure.jl and tests/test_kernel_regression.jl.
    bare_sol = solve(bare)
    @test maximum(abs, equilibrium_residuals(bare,
        [bare_sol.prices_raw; bare_sol.quantities])) < 1e-5
    # F1 is purely compositional and determinate on the open fixture too.
    f1 = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, 1.0;
        closure = :fixed, financing = PreferenceReallocation(fx.shift))
    f1_sol = solve(f1)
    @test maximum(abs, equilibrium_residuals(f1,
        [f1_sol.prices_raw; f1_sol.quantities])) < 1e-5
end

@testset "promoted closures: v3 data compatibility defaults" begin
    io = DataFrame("Sektoren" => ["a", "b", "c"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
    Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
    cs = [0.5, 0.3, 0.2]
    fs = [0.6, 0.5, 0.4]
    λ = [1.2, 1.0, 0.8]
    ls = λ .* fs
    data = Data(io, Ω, cs, fs, λ, ls, [0.4, 0.35, 0.3],
        [10.0, 8.0, 6.0], fs .* [10.0, 8.0, 6.0])
    @test data.gov_demand ≈ zeros(3)
    @test data.exo_demand ≈ zeros(3)
    @test data.exports_demand ≈ zeros(3)
    @test data.import_margin ≈ zeros(3)
    @test data.saving_rate ≈ 0.0 atol=1e-12
    @test data.household_baseline ≈ cs .* sum(ls)
end

@testset "promoted closures: Cobb-Douglas guard (verify_v3 §5)" begin
    fx = v3_fixture()
    shocks = _v3_shocks()
    fin = ExternalDebt(fx.g)
    ref = mobile_labor_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ, _V3_η;
        financing = fin)
    ref_sol = solve(ref)
    init = [ref_sol.prices_raw; ref_sol.quantities; ref_sol.wages_raw[1]; ref_sol.external_transfer]
    cd = mobile_labor_model(fx.data, shocks, 1.0, 1.0, 0.9, _V3_η;
        financing = fin)
    sol_cd = solve(cd; init = init)
    cd_near = mobile_labor_model(fx.data, shocks, 1.0, 1 - 1e-7, 0.9, _V3_η;
        financing = fin)
    sol_near = solve(cd_near; init = init)
    @test all(isfinite, sol_cd.prices_raw)
    @test all(isfinite, sol_near.prices_raw)
    @test abs(real_gdp(sol_cd) - real_gdp(sol_near)) < 1e-5
end

@testset "promoted closures: closure registry" begin
    @test closure_ids() ==
        [:BF, :ALPHA, :BETA, :GAMMA, :DELTA, :ZETA, :F1, :F2, :F3]
    for id in (:BF, :ALPHA, :BETA, :GAMMA, :DELTA, :ZETA)
        @test closure_axis(id) === :labor
    end
    for id in (:F1, :F2, :F3)
        @test closure_axis(id) === :financing
    end
    @test closure_constructor(:BF) === bf_model
    @test closure_constructor(:ALPHA) === alpha_model
    @test closure_constructor(:BETA) === beta_model
    @test closure_constructor(:GAMMA) === gamma_model
    @test closure_constructor(:DELTA) === delta_model
    @test closure_constructor(:ZETA) === nothing
    @test closure_constructor(:F1) === PreferenceReallocation
    @test closure_constructor(:F2) === TaxFinanced
    @test closure_constructor(:F3) === ExternalDebt
    @test_throws ArgumentError closure_constructor(:NOPE)
    @test_throws ArgumentError closure_axis(:NOPE)
    # The labour builders construct the advertised closures.
    fx = v3_fixture()
    shocks = _v3_shocks()
    @test labor_closure(bf_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
        0.5)) isa FlexibleWageClosure
    @test labor_closure(alpha_model(fx.data, shocks, _V3_θ, _V3_ϵ,
        _V3_σ)) isa FlexibleWageClosure
    @test labor_closure(beta_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
        _V3_η; eta_s = 0.5)) isa ElasticLaborClosure
    @test labor_closure(gamma_model(fx.data, shocks, _V3_θ, _V3_ϵ, _V3_σ,
        _V3_η)) isa FixedWageClosure
end
