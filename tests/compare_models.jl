#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using CSV
using CairoMakie
using DataFrames
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
cd_options_ls = CobbDouglas(cd_elasticities)

cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
cd_options = CES(cd_elasticities, model -> data.labor_share)

sol_cd = solve(Model(data, shocks, cd_options_ls))
sol_cd = solve(Model(data, shocks, cd_options))

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

## Calculting the shock amount proportional to GDP
gdp_effect_simple = 1 + shock_amount / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])

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
	s = copy(sol)
	len = 1000
	a = ones(len)
	arr = copy(el)
	for (idx, i) in enumerate(range(0.99, 0.015, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		try
			s = solve(model, init = vcat(s.prices, s.quantities))
			if labor_reallocation
				a[idx] = sol.gdp[1]
			else
				a[idx] = nominal ? s |> nominal_gdp : s |> real_gdp
			end
		catch
			a[idx] = NaN
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
		lines!(ax[i], [0.9, 0.015], 100 .* fill(gdp(sol_leontief, model_leontief), 2), label = "Elasticity between goods", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(gdp_effect_simple, 2), label = "Elasticity between goods", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(real_gdp(cd), 2), label = "Elasticity between goods", linestyle = :dash)
	end

	f[1, 2] = Legend(f, ax[1], labelsize = 25)

	f
end

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


p1 = plot_elasticities([a, b, c], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour slack")
p2 = plot_elasticities([d, e, f], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, without labour slack")
p3 = plot_elasticities([g, h, i], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, with labour reallocation")
save("plots/elastictiy_gradient_ls.png", p1)
save("plots/elastictiy_gradient_no_ls.png", p2)
save("plots/elastictiy_gradient_lr.png", p3)



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
lines!(ax, [0; 100], 100 .* fill(sol_cd |> real_gdp, 2), label = "Cobb Douglas", linestyle = :dash)
f[1, 2] = Legend(f, ax)
f

save("plots/labor_slack_gradient.png", f)


#============================================================================= 
Testing grounds
===============================================================================#
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
