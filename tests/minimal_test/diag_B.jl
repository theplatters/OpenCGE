# Milestone B verification: real_gdp is now the Tornqvist index. Confirm
# residuals still ~1e-16 (A intact) and inspect the real-GDP response.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
function maxres(model, sol)
    maximum(abs, equilibrium_residuals(sol))
end

noshock = Shocks(ones(71), ones(71), zeros(71))
base_model = mobile_labor_model(data, noshock, 0.5, 0.5, 0.9, 0.0)
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)
@printf("baseline real_gdp (Tornqvist) = %.6f   nominal_gdp = %.6f\n", base_gdp, nominal_gdp(base_sol))

println("eta     maxRes        realGDP      realGDPΔ%     nominalGDP    wage")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    shocks = standard_shock(data)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol); n = nominal_gdp(sol)
    @printf("%6.2f  %.2e  %10.6f  %+.6f%%  %10.6f  %.6f\n",
        η, mr, g, 100*(g/base_gdp - 1), n, sol.wages[1])
end
println("DONE")
