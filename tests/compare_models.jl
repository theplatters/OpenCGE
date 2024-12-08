#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using CSV
using CairoMakie
using DataFrames
using StatsBase
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
ces_options_ls = CES(ces_elasticities, model -> data.labor_share, true)

#cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
#cd_options = CES(cd_elasticities, x -> data.labor_share, false)
cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)

cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
cd_options = CES(cd_elasticities, full_labor_slack)

sol_cd_ls = solve(Model(data, shocks, cd_options))

ces_labor_realloc = CES(CESElasticities(0.2, 0.6, 0.9), x -> data.labor_share, true)
sol_realloc = solve(Model(data, shocks, ces_labor_realloc))
sol_realloc |> real_gdp

data.consumption_share |> sum

model = Model(data, shocks, ces_options)
sol = solve(model)
real_gdp(sol)

#=============================================================================
Simulating the  shock with a Leontief Model
===============================================================================#
leontief = Leontief()

model_leontief = Model(data, shocks, leontief)
sol_leontief = solve(model_leontief)

## Calculting the shock amount proportional to GDP
gdp_effect_simple = 1 + shock_amount / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])

#=============================================================================
Simulating the effect of different elasticites
===============================================================================#
struct ElasticityGradientSolution
	ϵ::Vector{Float64}
	θ::Vector{Float64}
	σ::Vector{Float64}
	mean_prices_ϵ::Vector{Float64}
	mean_prices_θ::Vector{Float64}
	mean_prices_σ::Vector{Float64}
	elasticities::CESElasticities
	labor_realloc::Bool
	nominal::Bool
end

function gradient(shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)
	s = copy(sol)
	len = 1000
	gdp = ones(len)
	mean_prices = ones(len)
	arr = copy(el)
	@inbounds for (idx, i) in enumerate(range(0.99, 0.015, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		try
			s = solve(model, init = vcat(s.prices, s.quantities))
			if labor_reallocation
				gdp[idx] = sol.gdp[1]
			else
				gdp[idx] = nominal ? s |> nominal_gdp : s |> real_gdp
				mean_prices[idx] = mean(s.prices, weights(s.quantities))
			end
		catch e
			@warn e
			gdp[idx] = NaN
			mean_prices[idx] = NaN
		end
		@show idx, arr, gdp[idx] mean_prices[idx]
	end
	return (gdp, mean_prices)
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

	(a, ap) = fetch(t1)
	(b, bp) = fetch(t2)
	(c, cp) = fetch(t3)
	elasticities = starting_elasticities
	return ElasticityGradientSolution(a, b, c, ap, bp, cp, CESElasticities(elasticities...), labor_reallocation, false)
end


function plot_elasticities(results; title = "Real GDP", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 720), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "GDP", ytickformat = "{:.2f}%", title = "0.9"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.5"),
		Axis(ga[1, 3], ytickformat = "{:.2f}%", title = "0.1")]

	supertitle = Label(f[0, :], title, fontsize = 40, tellwidth = false)
	linkaxes!(ax[1], ax[2], ax[3])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.σ), label = "Elasticity of consumption")
		lines!(ax[i], [0.9, 0.015], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Leontief model", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(gdp_effect_simple, 2), label = "Baseline Effect", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(real_gdp(cd), 2), label = "Cobb Douglas", linestyle = :dash)
	end

	f[1, 2] = Legend(f, ax[1], labelsize = 25)

	f
end

function plot_prices(results; title = "Real GDP", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 720), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "GDP", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]

	supertitle = Label(f[0, :], title, fontsize = 40, tellwidth = false)
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_σ), label = "Elasticity of consumption")
	end

	f[1, 2] = Legend(f, ax[1], labelsize = 25)

	f
end

a = elasticity_gradient(shocks, full_labor_slack, false, [0.99, 0.99, 0.99])
b = elasticity_gradient(shocks, full_labor_slack, false, [0.7, 0.7, 0.7])
c = elasticity_gradient(shocks, full_labor_slack, false, [0.2, 0.2, 0.2])
d = elasticity_gradient(shocks, full_labor_slack, false, [0.05, 0.05, 0.05])

e = elasticity_gradient(shocks, model -> data.labor_share, false, [0.99, 0.99, 0.99])
f = elasticity_gradient(shocks, model -> data.labor_share, false, [0.7, 0.7, 0.7])
g = elasticity_gradient(shocks, model -> data.labor_share, false, [0.2, 0.2, 0.2])
h = elasticity_gradient(shocks, model -> data.labor_share, false, [0.05, 0.05, 0.05])

## Labour reallocation
g = elasticity_gradient(shocks, model -> data.labor_share, true)
h = elasticity_gradient(shocks, model -> data.labor_share, true, [0.5, 0.5, 0.5])
i = elasticity_gradient(shocks, model -> data.labor_share, true, [0.1, 0.1, 0.1])


p1 = plot_elasticities([a, b, c], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour slack")
p2 = plot_elasticities([d, e, f], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, without labour slack")
p3 = plot_elasticities([g, h, i], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour reallocation")
save("plots/elastictiy_gradient_ls.png", p1)
save("plots/elastictiy_gradient_no_ls.png", p2)
save("plots/elastictiy_gradient_lr.png", p3)


p1 = plot_prices([a, b, c, d], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour slack")
p1 = plot_prices([e, f, g, h], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour slack")

#=============================================================================
Simulating labour slack effect
===============================================================================#

ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> full_labor_slack(model), false)
model = Model(data, shocks, ces)
sol = solve(model)
labour_slack_gradient = []
labour_slack_gradient_nominal = []
l(α, model) = (1 - α) * full_labor_slack(model) + α * model.data.labor_share
for α in range(0, 1, 100)
	labour_share(model) = l(α, model)
	ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> l(α, model))
	global model = Model(data, shocks, ces)
	global sol = solve(model, init = vcat(sol.prices, sol.quantities))
	push!(labour_slack_gradient, sol |> real_gdp)
	push!(labour_slack_gradient_nominal, sol |> nominal_gdp)
end

#============================================================================= 
Plotting the labour slack effect
===============================================================================#
f = Figure()
ax = Axis(f[1, 1], ytickformat = "{:.2f}%", ylabel = "GDP", xlabel = "Labour slack")
lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient, label = "Real GDP")
lines!(ax, [0; 100], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Leontief", linestyle = :dash)
lines!(ax, [0; 100], 100 .* fill(gdp_effect_simple, 2), label = "Baseline", linestyle = :dash)
lines!(ax, [0; 100], 100 .* fill(sol_cd_ls |> real_gdp, 2), label = "Cobb Douglas", linestyle = :dash)
f[1, 2] = Legend(f, ax)
f

save("plots/labor_slack_gradient.png", f)


#============================================================================= 
Testing grounds
===============================================================================#
cd_elasticities = CobbDouglasElasticities(data.factor_share, 1 .- data.factor_share)
ces_elasticities = CESElasticities(0.99, 0.99, 0.99)
options = CES(ces_elasticities, model -> full_labor_slack(model), false)
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
