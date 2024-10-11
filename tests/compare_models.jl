#============================================================================= 
Imports
===============================================================================#
using BeyondHulten
using DataFrames
using CSV
using Measures
using CairoMakie

#============================================================================= 
Loading in Data
===============================================================================#
data = Data("I-O_DE2019_formatiert.csv")

#============================================================================= 
Setting the shocks
===============================================================================#
shock_amount = 50000
demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [shock_amount]
calculate_investment!(shocks, data, investment, sector)


#============================================================================= 
Solving for commonly used elasticites and options
===============================================================================#
ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
ces_options = CES(ces_elasticities, x -> data.labor_share, false)
ces_options_ls = CES(ces_elasticities, data.labor_share, true)

cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
cd_options = CES(cd_elasticities, x -> data.labor_share, false)
cd_options_ls = CES(cd_elasticities, full_labor_slack, false)

sol_cd = solve(Model(data, shocks, cd_options))
sol_cd_ls = solve(Model(data, shocks, cd_options_ls))

ces_labor_realloc = CES(CESElasticities(0.2, 0.6, 0.9), x -> data.labor_share, true)
sol_realloc = solve(Model(data, shocks, ces_labor_realloc))
sol_realloc |> real_gdp



model = Model(data, shocks, ces_options)
sol = solve(model)
real_gdp(sol)

#============================================================================= 
Simulating the  shock with a Leontief Model
===============================================================================#
leontief = Leontief()

model_leontief = Model(data, shocks, leontief)
sol_leontief = solve(model_leontief)


#============================================================================= 
Simulating the effect of different elasticites
===============================================================================#
struct ElasticityGradientSolution
	ϵ::Vector{Float64}
	θ::Vector{Float64}
	σ::Vector{Float64}
	elasticities::CESElasticities
	labor_realloc::Bool
	nominal::Bool
end

function gradient(shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)

	len = 2000
	a = ones(len)
	arr = copy(el)
	for (idx, i) in enumerate(range(0.99, 0.001, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		if labor_reallocation
			a[idx] = sol.gdp[1]
		else
			a[idx] = nominal ? sol |> nominal_gdp : sol |> real_gdp
		end
		@info idx, arr, a[idx]
	end
	return a
end

function elasticity_gradient(shocks,
	labor_slack = full_labor_slack,
	labor_reallocation = false,
	starting_elasticities = [0.99, 0.99, 0.99],
	nominal = false)

	elasticities = CESElasticities(starting_elasticities...)
	ces = CES(elasticities, labor_slack, labor_reallocation)
	model = Model(data, shocks, ces)
	sol_original = solve(model)
	t1 = @task gradient(shocks, labor_slack, labor_reallocation, 1, sol_original, starting_elasticities, nominal)
	schedule(t1)
	t2 = @task gradient(shocks, labor_slack, labor_reallocation, 2, sol_original, starting_elasticities, nominal)
	schedule(t2)
	t3 = @task gradient(shocks, labor_slack, labor_reallocation, 3, sol_original, starting_elasticities, nominal)
	schedule(t3)

	a = fetch(t1)
	b = fetch(t2)
	c = fetch(t3)
	elasticities = starting_elasticities
	return ElasticityGradientSolution(a, b, c, CESElasticities(elasticities...), labor_reallocation, false)
end

## Calculting the shock amount proportional to GDP
gdp_effect_simple = 1 + shock_amount / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])

function plot_elasticities(results; title = "Real GDP", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 720))

	ga = f[1,1] = GridLayout()

	ax = [Axis(ga[1,1],  ylabel = "GDP", ytickformat = "{:.2f}%" ),
		Axis(ga[1,2], xlabel = "Elasticity", ytickformat = "{:.2f}%"),
		Axis(ga[1,3], ytickformat = "{:.2f}%")]

	linkaxes!(ax[1],ax[2],ax[3])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.025..0.9, 100 .* reverse(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.025..0.9, 100 .* reverse(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.025..0.9, 100 .* reverse(el.σ), label = "Elasticity of consumption")
		lines!(ax[i], [0.9, 0.025], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Elasticity between goods", linestyle = :dash)
		lines!(ax[i], [0.9, 0.025], 100 .* fill(gdp_effect_simple, 2), label = "Elasticity between goods", linestyle = :dash)
		lines!(ax[i], [0.9, 0.025], 100 .* fill(real_gdp(cd), 2), label = "Elasticity between goods", linestyle = :dash)
	end

	f[1,2] = Legend(f,ax[1], )

	f
end


demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [shock_amount]
calculate_investment!(shocks, data, investment, sector)

a = elasticity_gradient(shocks, full_labor_slack, false)
b = elasticity_gradient(shocks, full_labor_slack, false, [0.5, 0.5, 0.5])
c = elasticity_gradient(shocks, full_labor_slack, false, [0.1, 0.1, 0.1])


d = elasticity_gradient(shocks, model -> data.labor_share, false)
e = elasticity_gradient(shocks, model -> data.labor_share, false, [0.5, 0.5, 0.5])
f = elasticity_gradient(shocks, model -> data.labor_share, false, [0.1, 0.1, 0.1])

## Labour reallocation
g = elasticity_gradient(shocks, model -> data.labor_share, true)
h = elasticity_gradient(shocks, model -> data.labor_share, true, [0.5, 0.5, 0.5])
i = elasticity_gradient(shocks, model -> data.labor_share, true, [0.1, 0.1, 0.1])


p1 = plot_elasticities([a, b, c], cd = sol_cd_ls)
p2 = plot_elasticities([d, e, f], cd = sol_cd_ls)
p3 = plot_elasticities([g, h, i], cd = sol_cd_ls)
save("plots/makie_test.png",p1)


legend = plot([0 0 0 0 0 0], showaxis = false, grid = false, label = [" Elasticity between goods" " Elasticity between goods and labour" " Elasticity of consumption" " Baseline" " Leontief" " Cobb Douglas"])
legend_filler = plot([0 0 0 0], showaxis = false, grid = false, legend = :none, legendfontsize = 1000)
pall = Plots.plot(p1, p5, p9, legend_filler, legend, legend_filler,
	layout = @layout([grid(1, 3); [D{0.75w} E F{0.5h}]]),
	plot_title = "Effect of different elasticities on GDP, with labour slack",
	margin = 7mm,
	size = (1980, 1500),
	xtickfontsize = 14,
	ytickfontsize = 14,
	xguidefontsize = 14,
	yguidefontsize = 14,
	legendfontsize = 18,
	titlefontsize = 20,
	plot_titlefontsize = 24,
	legend_columns = 2)
Plots.savefig("plots/elasticity_gradient_50000_099.png")


pall = Plots.plot(p3, p7, p11, legend_filler, legend, legend_filler,
	layout = @layout([grid(1, 3); [D{0.75w} E F{0.5h}]]),
	plot_title = "Effect of different elasticities on GDP, without labour slack",
	margin = 7mm,
	size = (1980, 1500),
	xtickfontsize = 14,
	ytickfontsize = 14,
	xguidefontsize = 14,
	yguidefontsize = 14,
	legendfontsize = 18,
	titlefontsize = 20,
	plot_titlefontsize = 24,
	legend_columns = 2)
Plots.savefig("plots/elasticity_gradient_labour_slack.png")

pall = Plots.plot(p13, p15, p17, legend_filler, legend, legend_filler,
	layout = @layout([grid(1, 3); [D{0.75w} E F{0.5h}]]),
	plot_title = "Effect of different elasticities on GDP, without labour slack",
	margin = 7mm,
	size = (1980, 1500),
	xtickfontsize = 14,
	ytickfontsize = 14,
	xguidefontsize = 14,
	yguidefontsize = 14,
	legendfontsize = 18,
	titlefontsize = 20,
	plot_titlefontsize = 24,
	legend_columns = 2)
Plots.savefig("plots/elasticity_gradient_labour_realloc.png")

#============================================================================= 
Simulating labour slack effect
===============================================================================#
demand_shock = ones(71)
supply_shock = ones(71)
shocks = Shocks(demand_shock, supply_shock)
sector = ["Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"]
investment = [shock_amount]
calculate_investment!(shocks, data, investment, sector)
ces = CES(CESElasticities(0.01, 0.2, 0.9), model -> full_labor_slack(model), false)
model = Model(data, shocks, ces)
sol = solve(model)
labour_slack_gradient = []
labour_slack_gradient_nominal = []
l(α, model) = (1 - α) * full_labor_slack(model) + α * model.data.labor_share
for α in range(0, 1, 100)
	labour_share(model) = l(α, model)
	ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> l(α, model))
	model = Model(data, shocks, ces)
	sol = solve(model, init = vcat(sol.prices, sol.quantities))
	push!(labour_slack_gradient, sol |> real_gdp)
	push!(labour_slack_gradient_nominal, sol |> nominal_gdp)
end

f = Figure()
Axis(f[1, 1])
lines!(f[1, 1], range(100, 0, 100), 100 .* labour_slack_gradient, label = "Real GDP", color = :blue)
lines!(f[1, 1], [0; 100], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Real GDP", color = :red)
lines!(f[1, 1], [0; 100], 100 .* fill(gdp_effect_simple, 2), label = "Real GDP", color = :green)
lines!(f[1, 1], [0; 100], 100 .* fill(gdp(sol3_leontief, model_leontief), 2), label = "Real GDP", color = :red)

f
## Plotting the results
p1 = plot(range(100, 0, 100), 100 .* labour_slack_gradient, title = "Effect of different levels of labour slack on GDPs", label = "Effect on GDP", ylabel = "% of GDP", legend = :outerbottom)
plot!([0; 100], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model", linestyle = :dash, legendfontsize = 8)
plot!([0; 100], 100 .* fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP", linestyle = :dot)
plot!([0; 100], 100 .* fill(sol_cd |> real_gdp, 2), label = "Cobb Douglas GDP", linestyle = :dot)
plot!([0; 100], 100 .* fill(sol_cd |> real_gdp, 2), label = "Cobb Douglas GDP", linestyle = :dot)

p2 = plot(range(0, 1, 100), labour_slack_gradient_nominal, label = "Nominal GDP")
plot!([0; 1], fill(gdp(sol_leontief, model_leontief), 2), label = "Calculated GDP of the Leontief Model")
plot!([0; 1], fill(gdp_effect_simple, 2), label = "Proportion of shock to GDP")
plot!([0; 1], 1 .* fill(sol_cd |> real_gdp, 2), label = "Cobb Douglas GDP", linestyle = :dot)
plot(p1, p2, size = (1000, 500))

## Saving the plot

savefig("plots/labour_slack_gradient.png")

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
