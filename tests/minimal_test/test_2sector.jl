# Two-sector smoke diagnostic using the package implementation.
include("bootstrap.jl")

Ω = [0.6 0.4; 0.3 0.7]
consumption_share = [0.55, 0.45]
factor_share = [0.5, 0.6]
λ = transpose(I - diagm(1 .- factor_share) * Ω) \ consumption_share
labor_share = λ .* factor_share
data = Data(
    DataFrame(Sektoren=["sector 1", "sector 2"]),
    Ω,
    consumption_share,
    factor_share,
    λ,
    labor_share,
    ones(2),
    ones(2),
    factor_share,
)

supply = [1.1, 1.0]
shocks = Shocks(supply, ones(2), zeros(2))
println("η       real GDP    wage       max residual")
for η in [0.0, 0.5, 1.0]
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    @printf("%-7.2f %.6f    %.6f    %.2e\n", η, real_gdp(sol),
        sol.wages_raw[1], maximum(abs, equilibrium_residuals(sol)))
end
