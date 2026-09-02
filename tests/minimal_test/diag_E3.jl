# Test: does an AGGREGATE (uniform) autonomous demand shock engage the
# mobile-labor bridge (eta-dependence), whereas the sectoral shock does not?
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")

function maxres(model, sol)
    maximum(abs, equilibrium_residuals(sol))
end

# Baseline
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)

# Aggregate (uniform) autonomous demand shock: 0.5% across all sectors.
agg_shocks = Shocks(ones(71), ones(71), zeros(71), fill(0.005, 71), zeros(71))
println("=== AGGREGATE autonomous demand shock (0.5% uniform) ===")
println("eta     maxRes        realGDPΔ%     wage")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    model = mobile_labor_model(data, agg_shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f\n", 0.0, mr, 100*(g/base_gdp - 1), sol.wages[1])
end
# (loop above prints eta=0.0 always due to formatting bug; fix:)
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    model = mobile_labor_model(data, agg_shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol)
    @printf("eta=%6.2f  maxRes=%.2e  realGDPΔ=%+.6f%%  wage=%.6f\n", η, mr, 100*(g/base_gdp - 1), sol.wages[1])
end
println("DONE")
