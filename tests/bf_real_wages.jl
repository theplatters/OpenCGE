using BeyondHulten
using CSV
using GLMakie
using DataFrames
using StatsBase

data = Data("I-O_DE2019_formatiert.csv")
ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
ces_options = CES(ces_elasticities, x -> data.labor_share, false)


shock_range1 = range(1.0, 1.3, length = 50)

demand_shock = ones(71)
supply_shock = ones(71)
rgdp = ones(101)
ngdp = ones(101)
for (i,shock_amount) in enumerate(shock_range1)
	supply_shock[5] = shock_amount
	shocks = Shocks(demand_shock, supply_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	rgdp[i + 49] = sol.real_gdp[1]
	ngdp[i + 49] = sol.nominal_gdp[1]
end

shock_range2 = range(1.0, 0.7, length = 49)

demand_shock = ones(71)
supply_shock = ones(71)
gdp = ones(100)
for (i,shock_amount) in enumerate(shock_range2)
	supply_shock[5] = shock_amount
	shocks = Shocks(demand_shock, supply_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	rgdp[50 - i] = sol.real_gdp[1]
	ngdp[50 - i] = sol.nominal_gdp[1]
end

lines(1:99, ngdp[1:end-2], color = :blue)