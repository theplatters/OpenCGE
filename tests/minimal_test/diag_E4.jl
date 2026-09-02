include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)
println("baseline: wage=", base_sol.wages[1], "  sum(L_i)=", sum(sectoral_labor_demand(base_sol.prices, base_sol.quantities, base_sol.wages[1], base_model)), "  labor_bar=", base_model.options.labor_bar)

for m in [0.005, 0.03, 0.10]
    agg = Shocks(ones(71), ones(71), zeros(71), fill(m, 71), zeros(71))
    println("\n=== aggregate autonomous demand shock: $(100m)% uniform ===")
    for η in [1.0, 50.0]
        model = mobile_labor_model(data, agg, 0.5, 0.5, 0.9, η)
        sol = solve(model)
        L = sectoral_labor_demand(
            sol.prices_raw, sol.quantities, sol.wages_raw[1], model)
        @printf("eta=%6.2f  wage=%.6f  sum(L_i)=%.6f  realGDPΔ=%+.6f%%\n",
            η, sol.wages[1], sum(L), 100*(real_gdp(sol)/base_gdp - 1))
    end
end
println("DONE")
