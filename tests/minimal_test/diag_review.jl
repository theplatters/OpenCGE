# Regression diagnostic for the equilibrium defects identified by the second review.
# The package implementation is the source of truth for equations and residuals.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
shocks = standard_shock(data)
model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0)
sol = solve(model)

N = length(data.factor_share)
residuals = equilibrium_residuals(sol)
labor = sectoral_labor_demand(
    sol.prices_raw, sol.quantities, sol.wages_raw[1], model)
income = sol.wages_raw[1] * sum(labor)
expenditure = sum(sol.prices_raw .* sol.consumption)

println("=== Review regression diagnostic (71-sector mobile labor, η=0) ===")
println("household budget residual       = ", expenditure - income)
println("max |zero-profit residual|      = ", maximum(abs, residuals[1:N]))
println("max |goods-market residual|     = ", maximum(abs, residuals[N+1:2N-1]))
println("labor-market residual           = ", residuals[2N])
println("numeraire residual              = ", residuals[2N+1])
println("total max residual              = ", maximum(abs, residuals))
