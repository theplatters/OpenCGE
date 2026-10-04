using BeyondHulten, Test, LinearAlgebra
using Random

isdefined(Main, :v3_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# ADR-0019/ADR-0020 and `src/core/equilibrium.jl:534`: the mobile system is
# homogeneous of degree 1 in (p, w, F); the CPI = 1 row is the ONLY
# gauge-breaker. This pins that structure so the manuscript's "a nominal
# numeraire is harmless" claim is exercised by the suite rather than asserted
# only in prose (the DE-0011 probe is a weaker, demand-only re-scaling).
#
# What is actually tested (not a tautology):
#   (a) at a solution, re-gauging (p, w, F) -> k leaves the physical block
#       solved and moves ONLY the CPI pin (which reads k - 1, not 0);
#   (b) off the solution, the residual mapping is gauge-covariant: zero-profit
#       is degree 1 in (p, w), clearing and the real-wage labour row are degree
#       0, and the CPI pin is covariant (k*(r+1) - 1), exactly;
#   (c) the fixed-wage closure pins w = 1 INSTEAD of a CPI row, so scaling p
#       alone must break zero profit (the missing gauge freedom is the point).

@testset "numeraire invariance: mobile block is gauge-covariant" begin
    fx = v3_fixture()
    data = fx.data
    N = length(data.factor_share)
    shocks = _v3_shocks()
    Θ, Ε, Σ, Η = _V3_θ, _V3_ϵ, _V3_σ, _V3_η   # 0.5, 0.5, 0.9, 1.0

    # mobile η = 1 endpoint, with a live supply elasticity so the labour row is exercised
    mdl = mobile_labor_model(data, shocks, Θ, Ε, Σ, Η;
        financing = ExternalDebt(fx.g), eta_s = 0.5)
    sol = solve(mdl)
    X = [sol.prices_raw; sol.quantities; sol.wages_raw[1]; sol.external_transfer]
    @test maximum(abs, equilibrium_residuals(mdl, X)) < 1e-6

    # (a) re-gauging a solution keeps the physical block solved; only the pin moves
    for k in (0.5, 2.0)
        Xk = [k .* X[1:N]; X[N+1:2N]; k * X[2N+1]; k * X[2N+2]]
        rk = equilibrium_residuals(mdl, Xk)
        @test maximum(abs, rk[1:2N+1]) < 1e-6          # physical block stays 0
        @test rk[2N+2] ≈ k - 1 atol = 1e-12             # only the CPI pin detects k
    end

    # (b) off the solution, the residual mapping is gauge-covariant
    rng = MersenneTwister(7)
    Xp = max.(X .* (1 .+ 0.01 .* randn(rng, length(X))), 1e-3)
    r = equilibrium_residuals(mdl, Xp)
    for k in (0.5, 2.0)
        Xk = [k .* Xp[1:N]; Xp[N+1:2N]; k * Xp[2N+1]; k * Xp[2N+2]]
        rk = equilibrium_residuals(mdl, Xk)
        @test rk[1:N] ≈ k .* r[1:N]                        rtol = 1e-9 atol = 1e-9  # zero profit: degree 1
        @test rk[N+1:2N] ≈ r[N+1:2N]                       rtol = 1e-9 atol = 1e-9  # clearing: degree 0
        @test rk[2N+1] ≈ r[2N+1]                           rtol = 1e-9 atol = 1e-9  # real-wage labour: degree 0
        @test rk[2N+2] ≈ k * (r[2N+2] + 1) - 1             rtol = 1e-9 atol = 1e-9  # pin gauge-covariance
    end

    # (c) the fixed-wage closure pins w = 1 instead of a CPI row: scaling p alone
    #      must break zero profit (there is no gauge freedom to absorb it).
    mfix = Model(data, shocks,
        MobileLaborCES(MobileLaborCESElasticities(Θ, Ε, Σ, 1.0), 1.0, :fixed))
    sfix = solve(mfix)
    Xf = [sfix.prices_raw; sfix.quantities]
    Xfk = copy(Xf); Xfk[1:N] .*= 2.0
    rfk = equilibrium_residuals(mfix, Xfk)
    @test maximum(abs, rfk[1:N]) > 1e-3    # zero profit no longer holds
end
