using BeyondHulten, Test, LinearAlgebra
using TOML

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# ADR-0023 supply arm. The scalar-BETA identification cells (2N+2) carry the
# ln L / ln w signature used to recover the labour-supply elasticity eta_s.
# The old unit test merely asserted the model's own labour-supply equation
# (cbase2/review S3.4) and so was circular. This test instead pins three
# INDEPENDENT restrictions that the equation alone does not imply:
#   (i)   the wage is technology-determined and invariant to eta_s;
#   (ii)  the implied ln L / ln w recovers eta_s;
#   (iii) the implied elasticity is invariant to the shock-magnitude ladder
#         (m11/m12/m13) -- the genuine structural restriction;
#   (iv)  the ALPHA control (eta_s = 0) pins employment at Lbar even under the
#         supply-shock ladder (non-baseline real wage);
#   (v)   the re-solved numbers agree with the executed manifests to 1e-10
#         (golden regression, the same pattern as test_kernel_regression).
#
# Real-table blocks need the (gitignored) IO table; they skip when absent.

repo_root() = normpath(joinpath(@__DIR__, ".."))

@testset "supply-arm eta_s recovery (non-circular)" begin
    root = repo_root()
    data = isfile(joinpath(root, "data", "I-O_DE2019_formatiert.csv")) ?
        recalibrate_open(read_data("I-O_DE2019_formatiert.csv"; datadir = root);
            exo_scale = 1.0) : nothing
    if data === nothing
        @test_skip "needs data/I-O_DE2019_formatiert.csv"
    else
        N = length(data.factor_share)
        g = data.gov_demand
        fin = TaxFinanced(g)                 # F2, as in the executed manifests
        Θ, Ε, Σ, Η = 0.5, 0.5, 0.9, 1.0
        Lbar = sum(data.labor_share)

        wage_at_m12 = Float64[]  # wage per η_s, collected at the fixed A1 = 1.2 magnitude
        for η_s in (0.5, 1.0, 2.0, 5.0)
            implied = Float64[]
            for mag in (1.1, 1.2, 1.3)
                A = ones(N); A[1] = mag
                shocks = Shocks(A, ones(N), zeros(N))
                mdl = mobile_labor_model(data, shocks, Θ, Ε, Σ, Η;
                    financing = fin, eta_s = η_s)
                sol = solve(mdl)
                w = sol.wages_raw[1]
                L = sum(sectoral_labor_demand(sol.prices_raw, sol.quantities, w, mdl))
                # wage collected at A1 = 1.2, see below
                # implied elasticity from the equilibrium, not from the equation
                push!(implied, log(L / Lbar) / log(w))
                @test implied[end] ≈ η_s atol = 1e-6
                # golden against the executed manifest at A1 = 1.2
                if mag == 1.2
                    tag = η_s == 0.5 ? "e05" : η_s == 1.0 ? "e1" :
                          η_s == 2.0 ? "e2" : "e5"
                    mpath = joinpath(root, "runs",
                        "supply_etas_s1-BETA-F2-$(tag)-m12", "manifest.toml")
                    gm = TOML.parsefile(mpath)["metrics"]
                    @test w ≈ gm["wage"]       atol = 1e-10
                    @test L ≈ gm["employment"] atol = 1e-10
                    push!(wage_at_m12, w)  # identical across η_s here
                end
            end
            # the implied elasticity is invariant to the magnitude ladder
            @test maximum(abs, diff(implied)) < 1e-6
        end
        # the wage is technology-determined: invariant to η_s at a fixed magnitude
        # (the magnitude ladder itself moves the wage; that is a separate, real effect)
        @test maximum(abs, diff(wage_at_m12)) < 1e-10

        # ALPHA control: eta_s = 0 pins employment at Lbar (labour supply fixed).
        # Run through the SAME sector-1 supply-shock ladder as the recovery cells
        # above. At demand-only shocks (A = 1) the real wage stays at its
        # baseline, so any positive eta_s would also yield L = Lbar and the
        # assertion would be vacuous; under the ladder the real wage moves (see
        # the guard below) and a leaked positive supply elasticity would lift
        # employment by eta_s * log(w), failing the control.
        for mag in (1.1, 1.2, 1.3)
            A = ones(N); A[1] = mag
            mdlA = mobile_labor_model(data, Shocks(A, ones(N), zeros(N)),
                Θ, Ε, Σ, Η; financing = fin, eta_s = 0.0)
            solA = solve(mdlA)
            wA = solA.wages_raw[1]
            LA = sum(sectoral_labor_demand(solA.prices_raw, solA.quantities, wA, mdlA))
            @test LA ≈ Lbar atol = 1e-6
            # non-vacuity guard: the ladder moves the real wage by >= 2.6e-3 in
            # logs here, so a positive elasticity would be visible in LA
            @test abs(log(wA)) > 1e-3
        end
    end
end
