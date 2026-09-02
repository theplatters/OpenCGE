# Standalone package diagnostic for the mobile-labor model.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
N = length(data.factor_share)
baseline = solve(mobile_labor_model(
    data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 0.0))
shocks = standard_shock(data)

println("η        real GDP Δ%   wage       max residual")
for η in [0.0, 0.5, 1.0, 10.0, 100.0]
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    @printf("%-8.2f %+.6f%%    %.6f    %.2e\n", η,
        100 * (real_gdp(sol) / real_gdp(baseline) - 1),
        sol.wages_raw[1], maximum(abs, equilibrium_residuals(sol)))
end
