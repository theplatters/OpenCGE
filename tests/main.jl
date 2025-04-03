#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using GLMakie
using StatsBase
using ThreadsX
using ProgressMeter

const cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
const cd_options = CES(cd_elasticities, model -> model.data.labor_share)
const cd_options_ls = CES(cd_elasticities, BeyondHulten.full_labor_slack)
const ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
const ces_options = CES(ces_elasticities, model -> model.data.labor_share, false)
const ces_options_ls = CES(ces_elasticities, model -> model.data.labor_share, true)
const leontief = Leontief()


function plot_change_in_levels(data, impulses)
	shocks = impulse_shock(data, impulses)
	colors = Makie.wong_colors()
	sorted_shocks = sortperm(shocks.demand_shock, rev = true)  # Sort indices in descending order
	sorted_lambda = sortperm(data.λ, rev = true)  # Sort indices in descending order
	top5_indices = vcat(sorted_shocks[1:5], sorted_lambda[1:5]) #indices of sectors with biggest demand shock and with biggest overall size

	model = Model(data, shocks, ces_options)
	sol = solve(model)

	model_leontief = Model(data, shocks, Leontief())
	sol_leontief = solve(model_leontief)

	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], xlabel = "Sector",
		xticks = (1:length(top5_indices), data.io.Sektoren[top5_indices]),
		xticklabelrotation = -1 * pi / 4,
		ytickformat = "{:.2f}%")
	group = sort(repeat(1:4, length(top5_indices)))
	barplot!(ax,
		repeat(1:length(top5_indices), 4),
		[
			100 * (shocks.demand_shock[top5_indices] .- 1);
			100 * (shocks.demand_shock[top5_indices] .- 1);
			100 * (sol.consumption[top5_indices] ./ (data.consumption_share[top5_indices]) .- 1);
			100 * (sol.prices[top5_indices] .- 1)
		],
		dodge = group,
		color = colors[group])
	labels = ["Increase in state spending", "Change in consumption Leontief", "Change in consumption CGE", "Diviation of price from Numeraire CGE"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	axislegend(ax, elements, labels, position = :rt, labelsize = 20)
	save("plots/diff_consumption_imp.png", f)

	colors = Makie.wong_colors()
	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], xlabel = "Sector")
	shocks = impulse_shock(data, impulses)

	barplot!(ax,
		sol.prices .- 1,
		bar_labels = data.io.Sektoren[1:71], label_rotation = pi / 2, flip_labels_at = (0.0, 0.005))
	save("plots/diff_prices_imp.png", f)
end

function plot_change_price(data, impulse)
	shocks = impulse_shock(data, impulse)
	model = Model(data, shocks, ces_options)
	sol = solve(model)

	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1],
		xlabel = "Change in quantities",
		ylabel = "Change in prices",
		ytickformat = "{:.2f}%",
		xtickformat = "{:.2f}%")
	for (i, label) in enumerate(data.io.Sektoren[1:71])
		p = Point2f(100 * (sol.quantities[i] / data.λ[i] - 1), 100 * (sol.prices[i] - 1))
		scatter!(p, markersize = 20, color = :black)
		if (sol.prices[i] - 1 > 0.04)
			text!(p, text = label, color = :gray70, offset = (0, 20),
				align = (:center, :bottom))
		end
	end
	save("plots/scatter_prices_imp.png", f)
end

function diff_lambda(data, impulses)
	colors = Makie.wong_colors()
	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%")
	shocks = impulse_shock(data, impulses)

	model = Model(data, shocks, ces_options)
	sol = solve(model)
  sol_leontief = solve(Model(data, shocks, leontief))
	barplot!(ax,
    repeat(1:71, 2),
		[100 .* (sol.quantities ./ data.λ .- 1); 100 .* (sol_leontief.quantities[1:71] ./ data.λ .- 1)],
    stack = repeat(1:71,2),
    color = colors[sort(repeat(1:2, 71))])
	save("plots/diff_lambda_imp.png", f)
end


function simulate(shocks, data, title, gdp_effect_simple)
	a, b, c, d =
		ThreadsX.map([fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]) do elasticity
			BeyondHulten.elasticity_gradient(data, shocks, BeyondHulten.full_labor_slack, false, elasticity)
		end

	e, f, g, h =
		ThreadsX.map([fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]) do elasticity
			BeyondHulten.elasticity_gradient(data, shocks, model -> data.labor_share, false, elasticity)
		end

	sol_cd_ls = solve(Model(data, shocks, cd_options_ls))

	sol_cd = solve(Model(data, shocks, cd_options))

	model_leontief = Model(data, shocks, leontief)
	sol_leontief = solve(model_leontief)

	p1 = plot_real_gdp_gradient([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack ",
		cd = real_gdp(sol_cd_ls),
		leontief = gdp(sol_leontief, model_leontief),
		initial = gdp_effect_simple)

	p2 = plot_real_gdp_gradient([e, f, g, h],
		title = "Effect of different elasticities on GDP, without labour slack ",
		cd = real_gdp(sol_cd),
		leontief = real_gdp(sol_leontief),
		initial = gdp_effect_simple)


	p1_nominal = plot_nominal_gdp_gradient([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack ",
		cd = nominal_gdp(sol_cd_ls),
		leontief = nominal_gdp(sol_leontief),
		initial = gdp_effect_simple)

	p2_nominal = plot_nominal_gdp_gradient([e, f, g, h],
		title = "Effect of different elasticities on GDP, without labour slack ",
		cd = nominal_gdp(sol_cd),
		leontief = nominal_gdp(sol_leontief),
		initial = gdp_effect_simple)

	p1_wages = plot_wages([a, b, c, d],
		title = "Effect of different elasticities on wages with labour slack ",
		cd = real_gdp(sol_cd_ls))

	p2_wages = plot_wages([e, f, g, h],
		title = "Effect of different elasticities on wages without labour slack ",
		cd = real_gdp(sol_cd))

	p1_consumption = plot_consumption([a, b, c, d],
		title = "Effect of different elasticities on consumption without labour slack ",
		cd = real_gdp(sol_cd))
	p2_consumption = plot_consumption([e, f, g, h],
		title = "Effect of different elasticities on consumption without labour slack ",
		cd = real_gdp(sol_cd))


	save("plots/eg_$(title)_ls_rgdp.png", p1)
	save("plots/eg_$(title)_no_ls_rgdp.png", p2)

	save("plots/eg_$(title)_ls_ngdp.png", p1_nominal)
	save("plots/eg_$(title)_no_ls_ngdp.png", p2_nominal)

	save("plots/eg_$(title)_wages_ls.png", p1_wages)
	save("plots/eg_$(title)_wages_no_ls.png", p2_wages)

	save("plots/eg_$(title)_cons_ls.png", p1_consumption)
	save("plots/eg_$(title)_cons_no_ls.png", p2_consumption)
end

function plot_gradients(sectors, data)
	@showprogress Threads.@threads for sector in sectors
		shocks = Shocks(ones(71), ones(71))
		sector_number = findfirst(==(sector), data.io.Sektoren)
		shocks.demand_shock[sector_number] = 1.4
		gdp_effect_simple = 1 + 0.4 * data.io[sector_number, "Letzte Verwendung von Gütern zusammen"] / sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])
		simulate(shocks, data, sector, gdp_effect_simple)
	end
end

#=============================================================================
Simulating labour slack effect
===============================================================================#

function labor_slack_gradient(data, impulse)
	shocks = BeyondHulten.impulse_shock(data, impulse)
	gdp_effect_simple = 1 + sum(mean(col) for col in eachcol(impulse[:, 2:end-2] ./ sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"]')))
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	sol_cd_ls = solve(Model(data, shocks, cd_options))
	sol_leontief = solve(Model(data, shocks, leontief))
	labour_slack_gradient = Vector{Float64}()
	l(α, model) = (1 - α) * full_labor_slack(model) + α * model.data.labor_share
	for α in range(0, 1, 100)
		labour_share(model) = l(α, model)
		ces = CES(CESElasticities(0.01, 0.5, 0.9), model -> l(α, model))
		model = Model(data, shocks, ces)
		sol = solve(model, init = vcat(sol.prices, sol.quantities))
		push!(labour_slack_gradient, sol |> real_gdp)
	end

	f = Figure()
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%", ylabel = "GDP", xlabel = "Labour slack")
	lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient, label = "Real GDP")
	lines!(ax, [0; 100], 100 .* fill(real_gdp(sol_leontief), 2), label = "Leontief", linestyle = :dash)
	lines!(ax, [0; 100], 100 .* fill(gdp_effect_simple, 2), label = "Baseline", linestyle = :dash)
	lines!(ax, [0; 100], 100 .* fill(sol_cd_ls |> real_gdp, 2), label = "Cobb Douglas", linestyle = :dash)
	f[1, 2] = Legend(f, ax)
	f

	save("plots/labor_slack_gradient.png", f)
end



function (@main)(args)
	data = Data("I-O_DE2019_formatiert.csv")
	impulses = load_impulses("impulses.csv")
	sectors = [
		"Kohle",
		# "Chemische Erzeugnisse",
		"Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten",
		#"Fische, Fischerei- und Aquakulturerzeugnisse, DL",
		"Dienstleistungen v.Versicherungen u.Pensionskassen",
		"Nahrungs- u. Futtermittel, Getränke, Tabakerzeugn.",
		#"Sonstige Fahrzeuge",
	]


	plot_change_price(data, impulses)

	plot_gradients(sectors, data)
	shocks = impulse_shock(data, impulses)
	gdp_effect_simple = 1 + sum(mean(col) for col in eachcol(impulses[:, 2:end-2] ./ sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"]')))
	simulate(shocks, data, "impulse", gdp_effect_simple)
	plot_change_in_levels(data, impulses)
	diff_lambda(data, impulses)
	labor_slack_gradient(data, impulses)
end
