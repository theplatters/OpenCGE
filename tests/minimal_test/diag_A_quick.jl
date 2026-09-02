# Quick Milestone A check: baseline + a few eta values, foreground, file output.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
function full_residuals(model, sol)
    N = length(model.data.factor_share)
    residuals = equilibrium_residuals(sol)
    labor = sectoral_labor_demand(sol.prices_raw, sol.quantities,
        sol.wages_raw[1], model)
    budget = sum(sol.prices_raw .* sol.consumption) -
        sol.wages_raw[1] * sum(labor)
    return (zp=residuals[1:N], mc=residuals[N+1:2N-1],
        labor=residuals[2N], num=residuals[2N+1], budget=budget)
end
noshock = Shocks(ones(71), ones(71), zeros(71))
base_model = mobile_labor_model(data, noshock, 0.5, 0.5, 0.9, 0.0)
base_sol = solve(base_model); base_gdp = real_gdp(base_sol)
println("=== Quick Milestone A (corrected) — baseline + eta sweep ===")
println("eta      max|zp|      max|mc|      labor        num        budget      realGDPd%")
for η in [0.0, 0.5, 1.0]
    shocks = standard_shock(data)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model); r = full_residuals(model, sol); g = real_gdp(sol)
    @printf("%6.2f  %.2e  %.2e  %.2e  %.2e  %.2e  %+.3f%%\n",
        η, maximum(abs.(r.zp)), maximum(abs.(r.mc)), r.labor, r.num, r.budget, 100*(g/base_gdp - 1))
end
println("DONE")
