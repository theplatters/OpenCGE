# Deeper Milestone A check: precision + Tornqvist real-GDP index to see whether
# the mobile-labor bridge produces ANY eta-dependent real response.
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
base_sol = solve(base_model)
base_q = base_sol.quantities; base_p = base_sol.prices
base_gdp_laspeyres = real_gdp(base_sol)

println("eta     max|zp|      LaspeyresΔ%    Tornqvist(y)Δ%   wage")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    shocks = standard_shock(data)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_residuals(model, sol)
    g_las = 100*(real_gdp(sol)/base_gdp_laspeyres - 1)
    g_tor = 100*(tornqvist_quantity_index(sol.prices, sol.quantities, base_p, base_q) - 1)
    @printf("%6.2f  %.2e  %+.6f%%  %+.6f%%  %.6f\n",
        η, maximum(abs.(r.zp)), g_las, g_tor, sol.wages[1])
end
println("DONE")
