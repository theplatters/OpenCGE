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
rgdp38 = ones(100)
ngdp38 = ones(100)

rgdp35 = ones(100)
ngdp35 = ones(100)
for (i,shock_amount) in enumerate(shock_range1)
	supply_shock[35] = shock_amount
	shocks = Shocks(supply_shock,demand_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	rgdp35[i + 49] = sol.real_gdp[1]
	ngdp35[i + 49] = sol.nominal_gdp[1]
end

shock_range2 = range(1.0, 0.7, length = 49)

demand_shock = ones(71)
supply_shock = ones(71)
for (i,shock_amount) in enumerate(shock_range2)
	supply_shock[35] = shock_amount
	shocks = Shocks(supply_shock, demand_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
    if i != 1
        rgdp35[50 - i] = sol.real_gdp[1]
        ngdp35[50 - i] = sol.nominal_gdp[1]
    end
end

shock_range1 = range(1.0, 1.3, length = 50)
for (i,shock_amount) in enumerate(shock_range1)
	supply_shock[38] = shock_amount
	shocks = Shocks(supply_shock,demand_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	rgdp38[i + 49] = sol.real_gdp[1]
	ngdp38[i + 49] = sol.nominal_gdp[1]
end

shock_range2 = range(1.0, 0.7, length = 49)

demand_shock = ones(71)
supply_shock = ones(71)
for (i,shock_amount) in enumerate(shock_range2)
	supply_shock[38] = shock_amount
	shocks = Shocks(supply_shock, demand_shock)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
    if i != 1
        rgdp38[50 - i] = sol.real_gdp[1]
        ngdp38[50 - i] = sol.nominal_gdp[1]
    end
end

f = Figure()
ax = Axis(f[1,1], xlabel = "Shock", ylabel = "GDP")
lines!(ax,1:99, log.(ngdp35[1:end-1]), color = :blue)
lines!(ax,1:99, log.(rgdp35[1:end-1]), color = :red)

lines!(ax,1:99, log.(ngdp38[1:end-1]), color = :green)
lines!(ax,1:99, log.(rgdp38[1:end-1]), color = :yellow)
f