using OpenCGE, DataFrames, LinearAlgebra, Test

function _rw_data(d; kwargs...)
    overrides = Dict(kwargs)
    Data((get(overrides, f, getfield(d, f)) for f in fieldnames(Data))...)
end

# An open one-sector economy: y=1, L=1, c=0.8, investment=0.2 at baseline.
# With F=0 an extra g=0.1 requires y>1. Capacity/markup pricing then requires
# p/w>1, whereas a fixed real wage in one sector requires p/w=1. No nominal
# rescaling can remove that contradiction.
function _rw_one_sector()
    io = DataFrame("Sektoren" => ["a"],
        "Letzte Verwendung von Gütern zusammen" => [1.0])
    d = Data(io, ones(1, 1), [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0])
    _rw_data(d; saving_rate=0.2, household_baseline=[0.8], exo_demand=[0.2])
end

@testset "fixed real wage: constructors and demand-only employment adjustment" begin
    d = _rw_one_sector()
    sh = Shocks(ones(1), ones(1), zeros(1))
    for fin in (NoFinancing(), TaxFinanced([0.1]), ExternalDebt([0.1]))
        for m in (gamma_model(d, sh, 0.5, 0.5, 0.9, 1.0; financing=fin),
                  delta_model(d, sh; financing=fin))
            s = solve(m)
            @test labor_closure(m) isa FixedWageClosure
            @test s.wages_raw[1] / cpi(s) ≈ 1 atol=1e-10
            @test s.wages ≈ ones(1) atol=1e-10
            @test s.external_transfer == 0.0
            @test length(equilibrium_residuals(s)) == 3
            @test maximum(abs, equilibrium_residuals(s)) < 1e-9
            @test maximum(abs, market_clearing_residuals(m, [s.prices_raw; s.quantities])) < 1e-9
            for pr0 in (CapacityPricing(0.0), MarkupPricing(0.0))
                nested = solve(Model(d, sh, m.options, fin, pr0))
                @test nested.prices_raw == s.prices_raw
                @test nested.quantities == s.quantities
                @test nested.wages_raw == s.wages_raw
                @test nested.consumption == s.consumption
            end
            if fin isa ExternalDebt
                @test s.quantities ≈ [1.5] atol=1e-9
                @test sum(sectoral_labor_demand(s.prices_raw, s.quantities, 1.0, m)) ≈ 1.5 atol=1e-9
            end
        end
    end
    @test mobile_labor_model(d, sh, .5, .5, .9, 1.; closure=NominalWageClosure()).options.closure == :nominal
    @test labor_closure(MobileLaborCES(MobileLaborCESElasticities(.5, .5, .9, 1.), 1., :nominal)) isa NominalWageClosure
    @test_throws ArgumentError mobile_labor_model(d, sh, .5, .5, .9, 1.;
        closure=NominalWageClosure(), labor_bar=1.)
    m = gamma_model(d, sh, .5, .5, .9, 1.)
    @test labor_market_residual(FixedWageClosure(), m, 100., 3., 2.) == 0.5
    @test labor_market_residual(FixedWageClosure(), m, 100., 6., 4.) == 0.5
end

@testset "fixed real wage: incompatible pricing is rejected with F=0" begin
    d = _rw_one_sector()
    sh = Shocks(ones(1), ones(1), zeros(1))
    for pr in (CapacityPricing(0.5), MarkupPricing(0.5), MarkupPricing(0.5; s_pi=0.2))
        fin = ExternalDebt([0.1])
        for m in (gamma_model(d, sh, .5, .5, .9, 1.; financing=fin, pricing=pr),
                  delta_model(d, sh; financing=fin, pricing=pr))
            nominal = mobile_labor_model(d, sh, m.options.elasticities.θ,
                m.options.elasticities.ϵ, m.options.elasticities.σ, 1.;
                closure=NominalWageClosure(), financing=fin, pricing=pr)
            old = solve(nominal)
            @test abs(old.wages_raw[1] / cpi(old) - 1) > 1e-3
            @test old.external_transfer == 0.0
            x = [old.prices_raw; old.quantities]
            r = equilibrium_residuals(m, x)
            @test maximum(abs, r[1:2]) < 1e-9
            @test r[end] ≈ old.wages_raw[1] / cpi(old) - 1 atol=1e-12
            for init in (nothing, x)
                err = try
                    solve(m; init=init)
                catch e
                    e
                end
                @test err isa FixedRealWageFeasibilityError
                @test occursin("F = 0", sprint(showerror, err))
                @test occursin("goods residual", sprint(showerror, err))
            end
        end
    end
end

@testset "fixed real wage: uniform productivity shock has no compatible real wage" begin
    d = _rw_one_sector()
    # With fs=1, no pricing and A=1.1, zero profit requires p=w/1.1.
    # Thus w/CPI=1.1 at every nominal scale, independently of output.
    m = gamma_model(d, Shocks([1.1], ones(1), zeros(1)), .5, .5, .9, 1.)
    @test_throws FixedRealWageFeasibilityError solve(m)
    s = solve(mobile_labor_model(d, m.shocks, .5, .5, .9, 1.; closure=:nominal))
    @test s.wages_raw[1] / cpi(s) ≈ 1.1 atol=1e-9
    # Cobb–Douglas consumption needs its geometric CPI, not the degenerate
    # CES expression 1^Inf, which would silently remove the wage constraint.
    m_cd = gamma_model(d, m.shocks, .5, .5, 1., 1.)
    @test_throws FixedRealWageFeasibilityError solve(m_cd)
    s_cd = solve(mobile_labor_model(d, m.shocks, .5, .5, 1., 1.; closure=:nominal))
    @test cpi(s_cd) ≈ s_cd.prices_raw[1] atol=1e-12
    @test s_cd.wages_raw[1] / cpi(s_cd) ≈ 1.1 atol=1e-9
end

@testset "fixed real wage: feasible pricing with nonuniform prices" begin
    # Manufacture a compatible two-sector allocation with independent input
    # nests. Choose relative prices with CPI=1, output, and a pricing rule;
    # set productivity from unit cost, and investment from goods clearing.
    # The solver starts away from this allocation and must recover it.
    d0 = tiny_fixture()
    d0 = _rw_data(d0; saving_rate=0.4, household_baseline=[0.3, 0.3], exo_demand=zeros(2))
    y = [1.1, 0.9]
    for σ in (0.9, 1.0), ϵ in (0.5, 1.0), pr in (CapacityPricing(0.3), CapacityPricing(-0.3),
                                MarkupPricing(0.5), MarkupPricing(0.5; s_pi=0.2))
        p = [0.95, 1.05]
        p ./= σ == 1 ? sqrt(prod(p)) : sum(0.5 .* p .^ (1 - σ))^(1 / (1 - σ))
        μ = pr isa MarkupPricing ? 1 .+ pr.κ .* (y .- 1) : ones(2)
        unit_cost = ϵ == 1 ? sqrt.(p) : (0.5 .+ 0.5 .* sqrt.(p)) .^ 2
        A_eff = μ .* unit_cost ./ p
        A = pr isa CapacityPricing ? A_eff .* y .^ (pr.δ .* 0.5) : A_eff
        sh = Shocks(A, ones(2), zeros(2))
        trial = gamma_model(d0, sh, 0.5, ϵ, σ, 1.; pricing=pr)
        blocks = OpenCGE._mobile_market_demand(trial, p, y, 1.0)
        inv = y .- blocks.intermediary_demand .- blocks.total_final_demand
        @test all(>(0), inv)
        d = _rw_data(d0; exo_demand=inv)
        m = gamma_model(d, sh, 0.5, ϵ, σ, 1.; pricing=pr)
        s = solve(m; init=[1.03 .* p; 0.9 .* y])
        @test s.prices_raw ≈ p atol=1e-8
        @test s.quantities ≈ y atol=1e-8
        @test maximum(abs, s.prices_raw .- 1) > 0.04
        @test s.wages_raw[1] / cpi(s) ≈ 1 atol=1e-9
        @test s.numeraire ≈ cpi(s) atol=1e-12
        @test abs(sum(sectoral_labor_demand(s.prices_raw, s.quantities, 1., m)) - 1) > 1e-4
        @test s.external_transfer == 0.0
        @test maximum(abs, equilibrium_residuals(s)) < 1e-9
        @test abs(external_balance_canary(m, [s.prices_raw; s.quantities]).diff) < 1e-9
        adm = pricing_admissibility(m, s)
        @test adm.rank == adm.n == 4
        @test adm.converged == 4
        @test adm.spread < 1e-7
    end
end
