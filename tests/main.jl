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
const ces_options_ls = CES(ces_elasticities, BeyondHulten.full_labor_slack, false)
const leontief = Leontief()
const cd_options_ls_empirical = CES(cd_elasticities, BeyondHulten.empircial_labor_slack)
const ces_options_ls_empirical = CES(ces_elasticities, BeyondHulten.empircial_labor_slack)


function axis_change_in_level!(fig, data, impulses; options)
	shocks = impulse_shock(data, impulses)
	colors = Makie.wong_colors()
	sorted_shocks = sortperm(shocks.demand_shock, rev = true)  # Sort indices in descending order
	sorted_lambda = sortperm(data.λ, rev = true)  # Sort indices in descending order
	top5_indices = sorted_shocks[1:5]

	model = Model(data, shocks, options)
	sol = solve(model)

	model_leontief = Model(data, shocks, Leontief())
	sol_leontief = solve(model_leontief)

	# Error bar calculations
	model_low = Model(data, shocks, CES(CESElasticities(0.05, 0.05, 0.6), options.labor_slack, false))
	sol_low = solve(model_low)

	model_high = Model(data, shocks, CES(CESElasticities(0.99, 0.99, 0.99), options.labor_slack, false))
	sol_high = solve(model_high)

	ax = Axis(fig[1, 1], xlabel = "Sector",
		xticks = (1:length(top5_indices), data.io.Sektoren[top5_indices]),
		xticklabelrotation = -1 * pi / 4,
		ytickformat = "{:.2f}%")

	group = sort(repeat(1:4, length(top5_indices)))

	# Prepare the values for barplot
	consumption_values = 100 * (sol.consumption[top5_indices] ./ (data.consumption_share[top5_indices]) .- 1)
	price_values = 100 * (sol.prices[top5_indices] .- 1)

	# Prepare error ranges
	consumption_low = 100 * (sol_low.consumption[top5_indices] ./ (data.consumption_share[top5_indices]) .- 1)
	consumption_high = 100 * (sol_high.consumption[top5_indices] ./ (data.consumption_share[top5_indices]) .- 1)

	price_low = 100 * (sol_low.prices[top5_indices] .- 1)
	price_high = 100 * (sol_high.prices[top5_indices] .- 1)

	barplot!(ax,
		repeat(1:length(top5_indices), 4),
		[
			100 * (shocks.demand_shock[top5_indices] .- 1);
			100 * (sol_leontief.consumption[top5_indices] .- 1);
			consumption_values;
			price_values
		],
		dodge = group,
		color = colors[group])

	# Add error bars for the last two sets of bars (consumption and prices)
	for i in 1:length(top5_indices)
		# Position for consumption bars (3rd group)
		consumption_pos = i + (3 - 1) * 0.25 - 0.4
		# Position for price bars (4th group)
		price_pos = i + (4 - 1) * 0.25 - 0.45

		# Add error bars for consumption
		errorbars!(ax, [consumption_pos], [consumption_values[i]],
			[consumption_values[i] - consumption_low[i]],
			[consumption_high[i] - consumption_values[i]],
			whiskerwidth = 15, color = :black)

		# Add error bars for prices
		errorbars!(ax, [price_pos], [price_values[i]],
			[price_values[i] - price_low[i]],
			[price_high[i] - price_values[i]],
			whiskerwidth = 15, color = :black)
	end

	labels = ["Increase in state spending", "Change in consumption Leontief", "Change in consumption CGE", "Deviation of price from Numeraire CGE"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	axislegend(ax, elements, labels, position = :rt, labelsize = 20)
end

function get_color(data, shocks)
	colors = Makie.wong_colors()
	model = Model(data, shocks, ces_options)
	sol = solve(model)

	c = ones(Int, length(sol.prices))
	indices1 = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] < 1
	end

	indices2 = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] >= 1
	end
	c[indices1] .= 2
	c[indices2] .= 3
	c .+ 3
end
function axis_change_in_price!(fig, data, impulse; options)
	shocks = impulse_shock(data, impulse)
	model = Model(data, shocks, options)
	sol = solve(model)
	top5_indices = sortperm(abs.(sol.quantities ./ data.λ .- 1), rev = true)[1:5]

	colors = Makie.wong_colors()
	c = get_color(data, shocks)
	ax = Axis(fig[1, 2],
		xlabel = "Change in quantities",
		ylabel = "Change in prices",
		ytickformat = "{:.2f}%",
		xtickformat = "{:.2f}%")
	for (i, label) in enumerate(data.io.Sektoren[1:71])
		p = Point2f(100 * (sol.quantities[i] / data.λ[i] - 1), 100 * (sol.prices[i] - 1))
		scatter!(p, markersize = 20, color = colors[c[i]])
		if (i in top5_indices)
			text!(p, text = label, color = :gray70, offset = (0, 20),
				align = (:center, :bottom))
		end
	end
end

function panel(data, impulses; options = ces_options, name = "panel")
	fig = Figure(size = (1980, 1020))
	axis_change_in_level!(fig, data, impulses, options = options)
	axis_change_in_price!(fig, data, impulses, options = options)
	save("plots/$(name).png", fig)
end

function diff_lambda(data, impulses; options = ces_options, name = "diff_lambda_imp")
	colors = Makie.wong_colors()
	f = Figure(size = (1980, 1000))
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%")
	shocks = impulse_shock(data, impulses)

	model = Model(data, shocks, options)
	sol = solve(model)
	sol_leontief = solve(Model(data, shocks, leontief))
	grp = sort(repeat(1:2, 71))
	barplot!(ax,
		repeat(1:71, 2),
		[100 .* (sol.quantities ./ data.λ .- 1); 100 .* (sol_leontief.quantities[1:71] ./ data.λ .- 1)],
		dodge = grp,
		color = colors[grp])
	labels = ["CGE", "Leontief"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]

	axislegend(ax, elements, labels, position = :rt, labelsize = 28)
	save("plots/$(name).png", f)
end

function simulate(shocks, data, gdp_effect_simple; labor_slack_function = full_labor_slack, name = "eg")
	a, b, c, d =
		ThreadsX.map([fill(0.99, 3), fill(0.7, 3), fill(0.2, 3), fill(0.1, 3)]) do elasticity
			BeyondHulten.elasticity_gradient(data, shocks, labor_slack_function, false, elasticity)
		end


	sol_cd = solve(Model(data, shocks, CES(cd_elasticities, labor_slack_function)))
	model_leontief = Model(data, shocks, leontief)
	sol_leontief = solve(model_leontief)

	p1 = plot_real_gdp_gradient([a, b, c, d],
		title = "Effect of different elasticities on GDP with labour slack ",
		cd = real_gdp(sol_cd),
		leontief = real_gdp(sol_leontief),
		initial = gdp_effect_simple)

	@info gdp_effect_simple
	@info real_gdp(sol_cd)
	@info real_gdp(sol_leontief)
	save("plots/$(name).png", p1)
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
	sol_cd = solve(Model(data, shocks, cd_options))
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

	f = Figure(size = (1000, 800))
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%", ylabel = "GDP", xlabel = "Labour slack")
	lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient, label = "Real GDP")
	lines!(ax, [0, 100], 100 .* fill(real_gdp(sol_leontief), 2), label = "Leontief", linestyle = :dash)
	lines!(ax, [0, 100], 100 .* fill(real_gdp(sol_cd), 2), label = "Cobb-Douglas", linestyle = :dash)
	lines!(ax, [0, 100], 100 .* fill(gdp_effect_simple, 2), label = "Baseline", linestyle = :dot)
	axislegend(ax, position = :rb, labelsize = 28)
	save("plots/labor_slack_gradient.png", f)

	@info real_gdp(sol_cd)
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

	shocks = impulse_shock(data, impulses)
	gdp_effect_simple = 1 + sum(mean(col) for col in eachcol(impulses[:, 2:end-2] ./ sum(data.io[1:71, "Letzte Verwendung von Gütern zusammen"]')))

	panel(data, impulses)
	panel(data, impulses, options = ces_options_ls, name = "panel_ls")
	panel(data, impulses, options = ces_options_ls_empirical, name = "panel_ls_empirical")

	simulate(shocks, data, gdp_effect_simple, labor_slack_function = (model -> model.data.labor_share), name = "impulse")
	simulate(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.full_labor_slack, name = "impulse_ls")
	simulate(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.empircial_labor_slack, name = "impulse_ls_empirical")
	diff_lambda(data, impulses)
	diff_lambda(data, impulses, options = ces_options_ls, name = "diff_lambda_imp_lsl")
	diff_lambda(data, impulses, options = ces_options_ls_empirical, name = "diff_lambda_imp_ls_empirical")
	labor_slack_gradient(data, impulses)
end


function check_hunch(data, impulses)
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, ces_options)
	sol = solve(model)
	q = (data.Ω * sol.prices .^ (1 - sol.model.options.elasticities.θ)) .^ (1 / (1 - sol.model.options.elasticities.θ))

	indices = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] < 1
	end
	for (name, from_mean) in zip(data.io.Sektoren[indices], q[indices] .- mean(q)) |> collect
		println("Sector: ", name, " | q - mean(q): ", from_mean)
	end
	q[setdiff(1:71, indices)] |> mean
end