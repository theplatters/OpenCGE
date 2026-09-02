# Milestone A verification: after fixing the equilibrium system, confirm the
# solved allocation is a TRUE equilibrium (all residuals ~ 0) and re-run the
# eta-sweep to check the mobile-labor bridge no longer inverts.
include("bootstrap.jl")


data = Data("I-O_DE2019_formatiert.csv")

# Reconstruct the FULL corrected residual vector from a solved Solution.
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

# No-shock baseline for the Δ%-reference of real GDP.
noshock = Shocks(ones(71), ones(71), zeros(71))
base_model = mobile_labor_model(data, noshock, 0.5, 0.5, 0.9, 0.0)
base_sol   = solve(base_model)
base_gdp   = real_gdp(base_sol)

println("=== Milestone A verification (corrected system) ===")
println("η        max|zp|        max|mc|        labor          num          budget        realGDPΔ%")
ηs = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
for η in ηs
    shocks = standard_shock(data)
    model  = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    try
        sol = solve(model)
        r   = full_residuals(model, sol)
        g   = real_gdp(sol)
        @printf("%6.2f  %.2e  %.2e  %.2e  %.2e  %.2e  %+.3f%%\n",
                η, maximum(abs.(r.zp)), maximum(abs.(r.mc)),
                r.labor, r.num, r.budget, 100*(g/base_gdp - 1))
    catch e
        @printf("%6.2f  SOLVE FAILED: %s\n", η, e)
    end
end
