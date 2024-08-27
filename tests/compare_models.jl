using BeyondHulten
using Plots
data = Data("I-O_DE2019_formatiert.csv")
cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)
ces_elasticities = CESElasticities(0.001, 0.5, 0.9)


demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [50000]
calculate_investment!(shocks, data, investment, sector)


function elasticity_gradient(shocks)
	a = []
	b = []
	c = []

	a_nominal = []
	b_nominal = []
	c_nominal = []
	elasticities = CESELasticities(0.9, 0.5, 0.9)
	model = Model(data, shocks, cd_elasticities)
	sol = solve(model)
	for i in range(0.9, 0.1, 100)
		elasticities = CESELasticities(i, 0.5, 0.9)
		model = Model(data, shocks, elasticities)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		push!(a, (i, sol |> real_gdp))
		push!(a_nominal, (i, solve(model) |> nominal_gdp))
	end
	elasticities = CESELasticities(0.5, 0.9, 0.9)
	model = Model(data, shocks, cd_elasticities)
	sol = solve(model)
	for i in range(0.9, 0.1, 100)
		elasticities = CESELasticities(0.5, i, 0.9)
		model = Model(data, shocks, elasticities)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		push!(b, (i, sol |> real_gdp))
		push!(b_nominal, (i, sol |> nominal_gdp))
	end
	elasticities = CESELasticities(0.5, 0.5, 0.9)
	model = Model(data, shocks, cd_elasticities)
	sol = solve(model)
	for i in range(0.9, 0.1, 100)
		elasticities = CESELasticities(0.5, 0.5, i)
		model = Model(data, shocks, elasticities)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		push!(c, (i, sol |> real_gdp))
		push!(c_nominal, (i, sol |> nominal_gdp))
	end
	(a, b, c, a_nominal, b_nominal, c_nominal)
end


demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [25000]
calculate_investment!(shocks, data, investment, sector)

(a, b, c, d, e, f) = elasticity_gradient(shocks)

demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [50000]
calculate_investment!(shocks, data, investment, sector)

(a1, b1, c1, d1, e1, f1) = elasticity_gradient(shocks)

p1 = plot(first.(a), last.(a), label = "Goods", title = "Real GDP")
plot!(first.(b), last.(b), label = "Goods with Labour")
plot!(first.(c), last.(c), label = "Consumption")
p2 = plot(first.(d), last.(d), label = :none, title = "Nominal GDP")
plot!(first.(e), last.(e), label = :none)
plot!(first.(f), last.(f), label = :none)
p3 = plot(first.(a1), last.(a1), label = :none, title = "Real GDP")
plot!(first.(b1), last.(b1), label = :none)
plot!(first.(c1), last.(c1), label = :none)
p4 = plot(first.(d1), last.(d1), label = :none, title = "Nominal GDP")
plot!(first.(e1), last.(e1), label = :none)
plot!(first.(f1), last.(f1), label = :none)
plot(p1, p2, legend = true, size = (1000, 600))

savefig("plots/elasticity_gradient_25000.png")

model = Model(data, shocks, CESELasticities(0.01, 0.5, 0.9))
sol = solve(model)
labour_slack_gradient = []
labour_slack_gradient_nominal = []
l(α, model) = (1 - α) * full_demand_labor_allocation(model) + α * model.data.labor_share
for α in range(0, 1, 100)
	labour_share(model) = l(α, model)
	sol = solve(model, labor_reallocation = labour_share, init = vcat(sol.prices, sol.quantities))
	push!(labour_slack_gradient, sol |> real_gdp)
	push!(labour_slack_gradient_nominal, sol |> nominal_gdp)
end

p1 = plot(range(0,1,100),labour_slack_gradient, title = "Real GDP")
p2 = plot(range(0,1,100),labour_slack_gradient_nominal, title = "Nominal GDP")
plot(p1, p2, legend = false, size=(1000,500))

savefig("plots/labour_slack_gradient.png")
shocks.demand_shock[12] = 1.0

cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)
ces_elasticities = CESELasticities(0.001, 0.5, 0.9)
m2 = Model(data, shocks, cd_elasticities)
sol2 = solve(m2)

real_gdp(sol2)
nominal_gdp(sol2)	

using LinearAlgebra

elastities = LeontiefElasticiesLabor()
elastities = LeontiefElasticies()
demand_shock = Vector(impulses[1,3:end-1])
shocks = Shocks(demand_shock,demand_shock)

model = Model(data,shocks,elastities)
sol = solve(model)	

using CSV, DataFrames

impulses = CSV.read("data/impulses.csv",DataFrame)
