#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using CSV
using GLMakie
using DataFrames
using StatsBase
using ThreadsX
#=============================================================================
Loading in Data
===============================================================================#
data = Data("I-O_DE2019_formatiert.csv", read_unemployment = true)
impulses = load_impulses("impulses.csv")
#=============================================================================
Setting the shocks
===============================================================================#
shocks = standard_shock(data)
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

cd_options = CES(cd_elasticities, BeyondHulten.full_labor_slack_constrained)
sol_cd_ls2 = solve(Model(data, shocks, cd_options))
ces_labor_realloc = CES(CESElasticities(0.2, 0.6, 0.9), x -> data.labor_share, true)
sol_realloc = solve(Model(data, shocks, ces_labor_realloc))

#=============================================================================
Technology Shock
===============================================================================#
shocks = standard_tech_shock(data)

ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
ces_options = CES(ces_elasticities, x -> data.labor_share, false)
sol_ces_tech = solve(Model(data, shocks, ces_options))
#=============================================================================
Simulating the  shock with a Leontief Model
===============================================================================#
leontief = Leontief()

model_leontief = Model(data, shocks, leontief)
sol_leontief = solve(model_leontief)
shock_amount = 50_000
## Calculting the shock amount proportional to GDP

#=============================================================================
Simulating the effect of different elasticites
===============================================================================#

#shock 4 different sectors

sectors = [
	"Kohle",
	"Chemische Erzeugnisse",
	"Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten",
	"Fische, Fischerei- und Aquakulturerzeugnisse, DL",
	"Dienstleistungen v.Versicherungen u.Pensionskassen",
	"Nahrungs- u. Futtermittel, Getränke, Tabakerzeugn.",
	"Sonstige Fahrzeuge",
]

findfirst(==(sectors[1]),data.io.Sektoren)

shocks = impulse_shock(data, impulses)

effect = 1 .+ impulses[:, 2:end-2] ./ data.io[1:71, "Letzte Verwendung von Gütern zusammen"]'
for sector in sectors

	f = Figure(size = (2000, 2000))
	ax1 = Axis(f[1, 1], ylabel = "q - λ", xlabel = "Sector", title = "$(sector)")
	ax2 = Axis(f[2, 1], ylabel = "q - λ", xlabel = "sector -> shocked sector", title = "$(sector)")
	@info sector
	shocks = Shocks(ones(71), ones(71))
	sector_number = findfirst(==(sector), data.io.Sektoren)
	shocks.demand_shock[sector_number] = 1.8


	model = Model(data, shocks, ces_options)
	sol = solve(model)
	barplot!(ax1, sol.quantities - data.λ, bar_labels = data.io.Sektoren[1:71], label_rotation = pi / 2)
	scatter!(ax2, data.Ω[1:71, sector_number], sol.quantities - data.λ, markersize = 10, label = data.io.Sektoren[1:71], marker = :utriangle, color = 1:71, colormap = :plasma)
	scatter!(ax2, data.Ω[sector_number, 1:71], sol.quantities - data.λ, markersize = 10, label = data.io.Sektoren[1:71], marker = :dtriangle, color = 1:71, colormap = :plasma)
	save("plots/diff_lambda_$sector.png", f)
end


ces_options = CES(ces_elasticities, x -> data.labor_share, false)
begin
	colors = Makie.wong_colors()
	sorted_shocks = sortperm(shocks.demand_shock, rev=true)  # Sort indices in descending order
	sorted_lambda = sortperm(data.λ, rev = true)  # Sort indices in descending order
	top5_indices = vcat(sorted_indices[1:5],sorted_lambda[1:5])
	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], xlabel = "Sector", 
		xticks = (1:length(top5_indices), data.io.Sektoren[top5_indices]),
		xticklabelrotation = -1 * pi / 4)
	shocks = impulse_shock(data, impulses)

	model = Model(data, shocks, ces_options)
	sol = solve(model)

	model_leontief = Model(data, shocks, Leontief())
	sol_leontief = solve(model_leontief)
	#barplot!(ax, sol.quantities - data.λ, bar_labels=data.io.Sektoren[1:71], label_rotation=pi / 2, flip_labels_at=(0.0, 0.005))
	@info sum(sol.consumption .- data.consumption_share), real_gdp(sol)
	group = sort(repeat(1:4,length(top5_indices)))
	barplot!(ax,
		repeat(1:length(top5_indices), 4),
		[
			shocks.demand_shock[top5_indices] .- 1
			sol.quantities[top5_indices] ./ data.λ[top5_indices] .- 1;
			(data.λ[top5_indices] .+ sol_leontief.quantities_relative[top5_indices]) ./ data.λ[top5_indices] .- 1
			sol.prices[top5_indices] ./ mean(sol.prices,weights(sol.consumption)) .- 1
		],
		dodge = group,
		color = colors[group])
	labels = ["Increase in state spending", "Increase in quantities CGE", "Change in quantities Leontief", "Change in Price"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	title = "Groups"

	Legend(f[1, 2], elements, labels, title)
	save("plots/diff_consumption_imp.png", f)
end
f

begin
	colors = Makie.wong_colors()
	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], xlabel = "Sector")
	shocks = impulse_shock(data, impulses)

	model = Model(data, shocks, ces_options)
	sol = solve(model)
	#barplot!(ax, sol.quantities - data.λ, bar_labels=data.io.Sektoren[1:71], label_rotation=pi / 2, flip_labels_at=(0.0, 0.005))

	barplot!(ax,
		1 .- sol.prices ./ mean(sol.wages),
		bar_labels = data.io.Sektoren[1:71], label_rotation = pi / 2, flip_labels_at = (0.0, 0.005))
	save("plots/diff_prices_imp.png", f)
end



begin
	shocks = impulse_shock(data, impulses)
	@info shocks

	gdp_effect_simple = 1 + sum(mean(col) for col in eachcol(impulses[:, 2:end-2] ./ sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"]')))

	@info gdp_effect_simple
	a, b, c, d =
		fetch.([
			Threads.@spawn BeyondHulten.elasticity_gradient(data, shocks, full_labor_slack, false, elasticity) for elasticity in [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]])

	e, f, g, h =
		fetch.([
			Threads.@spawn BeyondHulten.elasticity_gradient(data, shocks, model -> model.data.labor_share, false, elasticity) for elasticity in [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]])


	cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
	cd_options = CES(cd_elasticities, full_labor_slack)
	sol_cd_ls = solve(Model(data, shocks, cd_options))
	leontief = Leontief()

	model_leontief = Model(data, shocks, leontief)
	sol_leontief = solve(model_leontief)

	#p1 = plot_prices([a, b, c, d], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, with labour slack " * sector)
	#p2 = plot_prices([e, f, g, h], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, without labour slack " * sector)
	p1 = plot_elasticities([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack ",
		cd = real_gdp(sol_cd_ls),
		leontief = gdp(sol_leontief, model_leontief),
		initial = gdp_effect_simple)

	p1_wages = plot_wages([a, b, c, d],
		title = "Effect of different elasticities on wages with labour slack ",
		cd = real_gdp(sol_cd_ls))

	cd_options = CES(cd_elasticities, model -> model.data.labor_share)
	sol_cd_ls = solve(Model(data, shocks, cd_options))

	p2 = plot_elasticities([e, f, g, h],
		title = "Effect of different elasticities on GDP, without labour slack ",
		cd = real_gdp(sol_cd_ls),
		leontief = gdp(sol_leontief, model_leontief),
		initial = gdp_effect_simple)

	p2_wages = plot_wages([e, f, g, h],
		title = "Effect of different elasticities on wages without labour slack ",
		cd = real_gdp(sol_cd_ls))

	save("plots/eg_imp_ls.png", p1)
	save("plots/eg_imp_no_ls.png", p2)

	save("plots/eg_imp_wages_ls.png", p1_wages)
	save("plots/eg_imp_wages_no_ls.png", p2_wages)
end



for sector in sectors
	@info sector
	shocks = Shocks(ones(71), ones(71))
	sector_number = findfirst(==(sector), data.io.Sektoren)
	shocks.demand_shock[sector_number] = 1.4
	gdp_effect_simple = 1 + 0.4 * data.io[sector_number, "Letzte Verwendung von Gütern zusammen"] / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])
	@info gdp_effect_simple
	a, b, c, d =
		ThreadsX.map(elasticity -> BeyondHulten.elasticity_gradient(data, shocks, full_labor_slack, false, elasticity), [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)])

	e, f, g, h =
		ThreadsX.map(elasticity -> BeyondHulten.elasticity_gradient(data, shocks, model -> data.labor_share, false, elasticity), [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)])


	cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
	cd_options = CES(cd_elasticities, full_labor_slack)
	sol_cd_ls = solve(Model(data, shocks, cd_options))
	leontief = Leontief()

	model_leontief = Model(data, shocks, leontief)
	sol_leontief = solve(model_leontief)

	#p1 = plot_prices([a, b, c, d], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, with labour slack " * sector)
	#p2 = plot_prices([e, f, g, h], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, without labour slack " * sector)
	p1 = plot_elasticities([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack " * sector,
		cd = real_gdp(sol_cd_ls),
		leontief = gdp(sol_leontief, model_leontief),
		initial = gdp_effect_simple)

	p1_wages = plot_wages([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack " * sector,
		cd = real_gdp(sol_cd_ls))

	cd_options = CES(cd_elasticities, model -> model.data.labor_share)
	sol_cd_ls = solve(Model(data, shocks, cd_options))

	p2 = plot_elasticities([e, f, g, h],
		title = "Effect of different elasticities on GDP, without labour slack " * sector,
		cd = real_gdp(sol_cd_ls),
		leontief = gdp(sol_leontief, model_leontief),
		initial = gdp_effect_simple)

	p2_wages = plot_wages([e, f, g, h],
		title = "Effect of different elasticities on GDP with labour slack " * sector,
		cd = real_gdp(sol_cd_ls))

	save("plots/elastictiy_gradient_ls_$sector.png", p1)
	save("plots/elastictiy_gradient_no_ls_$sector.png", p2)

	save("plots/elastictiy_gradient_wages_ls_$sector.png", p1_wages)
	save("plots/elastictiy_gradient_wages_no_ls_$sector.png", p2_wages)
end

for sector in sectors
	@info sector
	shocks = Shocks(ones(71), ones(71))

	sector_number = findfirst(==(sector), data.io.Sektoren)
	shocks.supply_shock[sector_number] = 1.1

	a, b, c, d =
		fetch.([
			Threads.@spawn BeyondHulten.elasticity_gradient(data, shocks, full_labor_slack, false, elasticity) for elasticity in [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]])

	e, f, g, h =
		fetch.([
			Threads.@spawn BeyondHulten.elasticity_gradient(data, shocks, model -> model.data.labor_share, false, elasticity) for elasticity in [fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]])


	cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
	cd_options = CES(cd_elasticities, full_labor_slack)
	sol_cd_ls = solve(Model(data, shocks, cd_options))
	leontief = Leontief()

	model_leontief = Model(data, shocks, leontief)
	sol_leontief = solve(model_leontief)
	p3 = plot_elasticities([a, b, c, d], cd = sol_cd_ls, title = "Effect of different elasticities on GDP with labour slack " * sector, cd = sol_cd_ls)
	p1 = plot_prices([a, b, c, d], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, with labour slack " * sector, cd = sol_cd_ls)

	cd_options = CES(cd_elasticities, data -> data.labor_share)
	sol_cd_ls = solve(Model(data, shocks, cd_options))
	p4 = plot_elasticities([e, f, g, h], cd = sol_cd_ls, title = "Effect of different elasticities on GDP, without labour slack " * sector, cd = sol_cd_ls)
	p2 = plot_prices([e, f, g, h], cd = sol_cd_ls, title = "Effect of different elasticities on mean prices, without labour slack " * sector, cd = sol_cd_ls)

	save("plots/eg_supply_ls_prices_$sector.png", p1)
	save("plots/eg_supply_no_ls_prices_$sector.png", p2)
	save("plots/eg_supply_ls_$sector.png", p3)
	save("plots/eg_supply_no_ls_$sector.png", p4)
end

#=============================================================================
concessions
===============================================================================#



for sector in sectors
	options = CES(CESElasticities(0.99, 0.99, 0.99), full_labor_slack, false)
	options_no_ls = CES(CESElasticities(0.99, 0.99, 0.99), model -> data.labor_share, false)
	demand_shock = standard_shock(data, sector)
	sol_demand = solve(Model(data, demand_shock, options))
	sol_demand_no_ls = solve(Model(data, demand_shock, options_no_ls))
	ts = standard_tech_shock(data, sector)
	@info ts.supply_shock[findfirst(==(sector), data.io.Sektoren)]
	sol_tech = solve(Model(data, ts, options))
	sol_tech_no_ls = solve(Model(data, ts, options_no_ls))



	concessions = DataFrame(
		sectors = data.io.Sektoren[1:71],
		price_deviation_demand = sol_demand.prices .- 1.0,
		price_deviation_tech = sol_tech.prices .- 1.0,
		price_deviation_demand_no_ls = sol_demand_no_ls.prices .- 1.0,
		price_deviation_tech_no_ls = sol_tech_no_ls.prices .- 1.0,
		from = Vector(data.io[1:71, sector]),
		to = Vector(data.io[findfirst(==(sector), data.io.Sektoren), 2:72]))

	CSV.write("data/concessions_$(sector).csv", concessions)
end
fig = Figure()
ax1 = Axis(fig[1, 1], ylabel = "Price Derivation", xlabel = "Sectors")
ax2 = Axis(fig[2, 1], ylabel = "Concessions", xlabel = "Sectors")

barplot!(ax1, concessions_tech.price_derivation, bar_labels = concessions_tech.sectors)
barplot!(ax2, concessions_tech.from, bar_labels = concessions_tech.sectors)


fig
#=============================================================================
Simulating labour slack effect
===============================================================================#


shocks = BeyondHulten.standard_shock(data)
ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> full_labor_slack(model), false)
model = Model(data, shocks, ces)
sol = solve(model)
sector_number = findfirst(==("Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"), data.io.Sektoren)
gdp_effect_simple = 1 + 0.4 * data.io[sector_number, "Letzte Verwendung von Gütern zusammen"] / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])
sol_cd_ls = solve(Model(data, shocks, cd_options))
labour_slack_gradient = Vector{Float64}()
labour_slack_gradient_prices = Vector{Float64}()
labour_slack_gradient_prices_weighted = Vector{Float64}()
labour_slack_gradient_nominal = Vector{Float64}()
l(α, model) = (1 - α) * full_labor_slack(model) + α * model.data.labor_share
for α in range(0, 1, 100)
	labour_share(model) = l(α, model)
	ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> l(α, model))
	global model = Model(data, shocks, ces)
	global sol = solve(model, init = vcat(sol.prices, sol.quantities))
	push!(labour_slack_gradient, sol |> real_gdp)
	push!(labour_slack_gradient_prices_weighted, mean(sol.prices, weights(sol.quantities)))
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

f = Figure()
ax = Axis(f[1, 1], ytickformat = "{:.2f}%", xlabel = "Labour slack")
lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient_prices, label = "Mean Prices")
lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient_prices_weighted, label = "Weighted Mean Prices")
f[1, 2] = Legend(f, ax)
f

save("plots/labor_slack_gradient_prices.png", f)


f = Figure()
ax = Axis(f[1, 1], ytickformat = "{:.2f}%", xlabel = "Labour slack")
lines!(ax, range(100, 0, 100), 100 .* (labour_slack_gradient - labour_slack_gradient_nominal), label = "Mean Prices")
f[1, 2] = Legend(f, ax)
f

#=============================================================================
Different Investments
===============================================================================#


using GLMakie

fig = Figure()

ax = Axis(fig[1, 1])

sl_x = Slider(fig[2, 1], range = 0.7:0.1:1.8, startvalue = 1.4)
sl_y = Slider(fig[1, 2], range = 1:71, startvalue = 1, horizontal = false)

ces_elasticities = CESElasticities(0.01, 0.5, 0.9)
ces_options = CES(ces_elasticities, full_labor_slack, false)
p = lift(sl_x.value, sl_y.value) do x, y
	@info data.io.Sektoren[y]
	demand_shock = ones(71)
	supply_shock = ones(71)
	demand_shock[y] = x
	shocks = Shocks(demand_shock, supply_shock)

	m = Model(data, shocks, ces_options)
	sol = solve(m)
	@info sol |> real_gdp, sol |> nominal_gdp
	Point2f(sol |> real_gdp, sol |> nominal_gdp)
end

scatter!(p, color = :red, markersize = 20)

limits!(ax, 0.9, 1.1, 0.9, 1.1)

fig
