using BeyondHulten
using Plots
using DataFrames
using CSV
using Measures

data = Data("I-O_DE2019_formatiert.csv")
cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)
ces_elasticities = CESElasticities(0.001, 0.5, 0.9)

ces_options = CES(ces_elasticities, x -> data.labor_share, true)
cd_options = CES(CESElasticities(0.99,0.99,0.99), x -> data.labor_share, true)
## setting the shock
shock_amount = 50000



demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [50000]
calculate_investment!(shocks, data, investment, sector)
sol_cd = solve(Model(data, shocks, cd_options))

sum(sol_cd.value_added_relative)

model = Model(data, shocks, ces_options)
sol = solve(model)
real_gdp(sol)
## Simulating the  shock with a Leontief Model
leontief = Leontief()

model_leontief = Model(data, shocks, leontief)
sol_leontief = solve(model_leontief)


function gradient(shocks, labor_slack, labor_reallocation, elasticity, sol, el)

	len = 200
	a = ones(len)
	a_nominal = ones(len)
	arr = [0.9, 0.9, 0.9]
	for (idx, i) in enumerate(range(0.90, 0.025, len))
		arr .= el
		arr[elasticity] = i
		@info idx, arr
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		a[idx] = sol |> real_gdp
		a_nominal[idx] = sol |> nominal_gdp
	end
	(a, a_nominal)
end

function elasticity_gradient(shocks,
	labor_slack = BeyondHulten.full_demand_labor_allocation,
	labor_reallocation = false,
	starting_elasticities = [0.99, 0.99, 0.99])

	elasticities = CESElasticities(starting_elasticities...)
	ces = CES(elasticities, labor_slack, labor_reallocation)
	model = Model(data, shocks, ces)
	sol_original = solve(model)
	t1 = @task gradient(shocks, labor_slack, labor_reallocation, 1, sol_original, starting_elasticities)
	schedule(t1)
	t2 = @task gradient(shocks, labor_slack, labor_reallocation, 2, sol_original, starting_elasticities)
	schedule(t2)
	t3 = @task gradient(shocks, labor_slack, labor_reallocation, 3, sol_original, starting_elasticities)
	schedule(t3)

	(a, a_nominal) = fetch(t1)
	(b, b_nominal) = fetch(t2)
	(c, c_nominal) = fetch(t3)
	elasticities = starting_elasticities
	sol = sol_original
	(a, b, c, a_nominal, b_nominal, c_nominal)
end

## Calculting the shock amount proportional to GDP
gdp_effect_simple = 1 + shock_amount / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])

function plot_elasticities(a, b, c, d, e, f; ran = range(0.9, 0.025, length(a)), title = "Real GDP")

	plot_gdp = plot(ran, 100 .* a, title = title, ylabel="%", ylims=(97,103), xlabel="Elasticity", xlims=(0,1), legend=:none)
	plot!(ran, 100 .* b, label = "Elasticity between labour and goods")
	plot!(ran,100 .*  c, label = "Elasticity of consumption")
	plot!([0.9, 0.025],100 .*  fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP", linestyle=:dot)
	plot!([0.9, 0.025],100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model", linestyle=:dash)
	plot_gdp_nominal = plot(ran, 100 .* d, label = "", title = title, ylabel="%", ylims=(97,103), xlabel="Elasticity", xlims=(0,1))
	plot!(ran, 100 .* e, label = "Elasticity between labour and goods")
	plot!(ran,100 .*  f, label = "Elasticity of consumption")
	plot!([0.9, 0.025],100 .*  fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP", linestyle=:dot)
	plot!([0.9, 0.025],100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model", linestyle=:dash)
	return plot_gdp, plot_gdp_nominal
end

demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [shock_amount]
calculate_investment!(shocks, data, investment, sector)

(a1, b1, c1, d1, e1, f1) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false)
(a2, b2, c2, d2, e2, f2) = elasticity_gradient(shocks, model -> data.labor_share, false)
(a3, b3, c3, d3, e3, f3) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false, [0.5,0.5,0.5])
(a4, b4, c4, d4, e4, f4) = elasticity_gradient(shocks, model -> data.labor_share, false, [0.5, 0.5, 0.5])
(a5, b5, c5, d5, e5, f5) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false,[0.1, 0.1, 0.1])
(a6, b6, c6, d6, e6, f6) = elasticity_gradient(shocks, model -> data.labor_share, false, [0.1, 0.1, 0.1])

p1, p2 = plot_elasticities(a1, b1, c1, d1, e1, f1, title = "0.9")
p3, p4 = plot_elasticities(a2, b2, c2, d2, e2, f2, title="0.9")

p5, p6 = plot_elasticities(a3, b3, c3, d3, e3, f3, title = "0.5")
p7, p8 = plot_elasticities(a4, b4, c4, d4, e4, f4, title = "0.5")

p9, p10 = plot_elasticities(a5, b5, c5, d5, e5, f5, title="0.1")
p11, p12 = plot_elasticities(a6, b6, c6, d6, e6, f6, title="0.1")
legend = plot([0 0 0  0 0], showaxis = false, grid = false, label = [" Elasticity between goods" " Elasticity between goods and labour" " Elasticity of consumption" " Baseline" " Leontief"])
legend_filler = plot([0 0 0 0], showaxis = false, grid = false, legend=:none, legendfontsize=1000) 
pall = Plots.plot(p1, p5, p9, legend_filler, legend,legend_filler, 
	layout = @layout([grid(1,3); [D{0.75w} E F{0.5h}]]), 
	plot_title="Effect of different elasticities on GDP, with labour slack", 
	margin=7mm, 
	size=(1980,1500), 
	xtickfontsize=14,
	ytickfontsize=14,
	xguidefontsize=14,
	yguidefontsize=14,
	legendfontsize=18,
	titlefontsize=20,
	plot_titlefontsize=24,
	legend_columns=2)
Plots.savefig("plots/elasticity_gradient_50000_099.png")


pall = Plots.plot(p3, p7, p11, legend_filler, legend,legend_filler, 
	layout = @layout([grid(1,3); [D{0.75w} E F{0.5h}]]), 
	plot_title="Effect of different elasticities on GDP, without labour slack", 
	margin=7mm, 
	size=(1980,1500), 
	xtickfontsize=14,
	ytickfontsize=14,
	xguidefontsize=14,
	yguidefontsize=14,
	legendfontsize=18,
	titlefontsize=20,
	plot_titlefontsize=24,
	legend_columns=2)
Plots.savefig("plots/elasticity_gradient_labour_slack.png")



## Simulating labour slack effect
demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [shock_amount]
calculate_investment!(shocks, data, investment, sector)
ces = CES(CESElasticities(0.01, 0.2, 0.9), model -> full_demand_labor_allocation(model), false)
model = Model(data, shocks, ces)
sol = solve(model)
labour_slack_gradient = []
labour_slack_gradient_nominal = []
l(α, model) = (1 - α) * full_demand_labor_allocation(model) + α * model.data.labor_share
for α in range(0, 1, 100)
	labour_share(model) = l(α, model)
	ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> l(α, model))
	model = Model(data, shocks, ces)
	sol = solve(model, init = vcat(sol.prices, sol.quantities))
	push!(labour_slack_gradient, sol |> real_gdp)
	push!(labour_slack_gradient_nominal, sol |> nominal_gdp)
end


## Plotting the results
p1 = plot(range(100,0, 100),100 .* labour_slack_gradient, title = "Effect of different levels of labour slack on GDPs", label = "Effect on GDP", ylabel="% of GDP")
plot!([0; 100], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model", linestyle=:dash, legendfontsize=8)
plot!([0; 100], 100 .* fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP", linestyle=:dot)

p2 = plot(range(0, 1, 100), labour_slack_gradient_nominal, title = "Nominal GDP")
plot!([0; 1], fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model")
plot!([0; 1], fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP")
plot(p1, p2, size = (1000, 500))

## Saving the plot

savefig("plots/labour_slack_gradient.png")
shocks.demand_shock[12] = 1.0

### Testing area

cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)
ces_elasticities = CESElasticities(0.99, 0.99, 0.99)
options = CES(ces_elasticities, model -> full_demand_labor_allocation(model), false)
m2 = Model(data, shocks, options)
sol2 = solve(m2)
sol2 |> real_gdp

real_gdp(sol2)
nominal_gdp(sol2)



using LinearAlgebra

impulses = CSV.read("data/impulses.csv", DataFrame) ./ 10e5
elastities = LeontiefElasticiesLabor()
elastities = LeontiefElasticies()
demand_shock = Vector(impulses[1, 3:end-1])
maximum(eachcol(impulses))
calculate_investment!(shocks, data, Vector(impulses[1, 3:end-1]), data.io.Sektoren[1:71])
shocks = Shocks(demand_shock, demand_shock)

shocks.demand_shock

impulses[1, "Machinery"]
names(impulses)

model = Model(data, shocks, ces_elasticities)
sol = solve(model)


data.io.Sektoren
## Supply shockdemand_shock = ones(71)
supply_shock = ones(71)
supply_shock[4] = 0.8
demand_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)


(a1, b1, c1, d1, e1, f1) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false)
(a2, b2, c2, d2, e2, f2) = elasticity_gradient(shocks, model -> data.labor_share, false)
(a3, b3, c3, d3, e3, f3) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false, [0.5,0.5,0.5])
(a4, b4, c4, d4, e4, f4) = elasticity_gradient(shocks, model -> data.labor_share, false, [0.5, 0.5, 0.5])
(a5, b5, c5, d5, e5, f5) = elasticity_gradient(shocks, BeyondHulten.full_demand_labor_allocation, false,[0.1, 0.1, 0.1])
(a6, b6, c6, d6, e6, f6) = elasticity_gradient(shocks, model -> data.labor_share, false, [0.1, 0.1, 0.1])

p1, p2 = plot_elasticities(a1, b1, c1, d1, e1, f1, title = "0.9")

p3, p4 = plot_elasticities(a2, b2, c2, d2, e2, f2, title="0.9")

p5, p6 = plot_elasticities(a3, b3, c3, d3, e3, f3, title = "0.5")
p7, p8 = plot_elasticities(a4, b4, c4, d4, e4, f4, title = "0.5")

p9, p10 = plot_elasticities(a5, b5, c5, d5, e5, f5, title="0.1")
p11, p12 = plot_elasticities(a6, b6, c6, d6, e6, f6, title="0.1")
legend = plot([0 0 0  0 0], showaxis = false, grid = false, label = [" Elasticity between goods" " Elasticity between goods and labour" " Elasticity of consumption" " Baseline" " Leontief"])
legend_filler = plot([0 0 0 0], showaxis = false, grid = false, legend=:none, legendfontsize=1000) 
pall = Plots.plot(p1, p5, p9, legend_filler, legend,legend_filler, 
	layout = @layout([grid(1,3); [D{0.75w} E F{0.5h}]]), 
	plot_title="Effect of different elasticities on GDP, with labour slack", 
	margin=7mm, 
	size=(1980,1500), 
	xtickfontsize=14,
	ytickfontsize=14,
	xguidefontsize=14,
	yguidefontsize=14,
	legendfontsize=18,
	titlefontsize=20,
	plot_titlefontsize=24,
	legend_columns=2)