# Diagnostic for the package's canonical intersectoral-reallocation semantics.
# η=0 is the baseline allocation, η=1 is the cost-minimizing allocation, and
# values above one extrapolate the same geometric rule.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
N = length(data.factor_share)
base_shocks = Shocks(ones(N), ones(N), zeros(N))
base_model = mobile_labor_model(data, base_shocks, 0.5, 0.5, 0.9, 0.0)
base_gdp = real_gdp(solve(base_model))

shocks = standard_shock(data)
println("Canonical η sweep (package equilibrium):")
println("  η        Real GDP Δ%   Total Labor   Wage       max residual")
println("  " * "-"^68)
for η in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    labor = sum(sectoral_labor_demand(
        sol.prices_raw, sol.quantities, sol.wages_raw[1], model))
    @printf("  %-8.2f %+.6f%%    %.6f      %.6f   %.2e\n",
        η, 100 * (real_gdp(sol) / base_gdp - 1), labor,
        sol.wages_raw[1], maximum(abs, equilibrium_residuals(sol)))
end
