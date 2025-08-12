
function axis_change_in_level!(fig, data, impulses; options)
	shocks = impulse_shock(data, impulses)
	colors = Makie.wong_colors()
	sorted_shocks = sortperm(shocks.demand_shock, rev = true)  # Sort indices in descending order
	sorted_lambda = sortperm(data.λ, rev = true)  # Sort indices in descending order
	top5_indices = sorted_shocks[1:7]

	model = Model(data, shocks, options)
	sol = solve(model)

	model_leontief = Model(data, shocks, Leontief())
	sol_leontief = solve(model_leontief)

	# Error bar calculations
	model_low = Model(data, shocks, CES(CESElasticities(0.001, 0.05, 0.6), options.labor_slack, false))
	sol_low = solve(model_low)

	model_high = Model(data, shocks, CES(CESElasticities(0.99, 0.99, 0.99), options.labor_slack, false))
	sol_high = solve(model_high)

	ax = Axis(fig[1, 1], xlabel = "Sector",
		xticks = (1:length(top5_indices), data.io.Sektoren[top5_indices]),
		xticklabelrotation = -1 * pi / 4,
		ytickformat = "{:.2f}%",
		ylabelsize = 24,
		xlabelsize = 24,
		xticklabelsize = 20,
		yticklabelsize = 20)

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
	axislegend(ax, elements, labels, position = :rt, labelsize = 24)
end

function get_color(data, shocks)
	colors = Makie.wong_colors()
    ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
    ces_options = CES(ces_elasticities, model -> model.data.labor_share, false)
	model = Model(data, shocks, ces_options)
	sol = solve(model)

	c = ones(Int, length(sol.prices))
	indices1 = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] < 1
	end

	top7_indices = sortperm(shocks.demand_shock, rev = true)[1:7]
	indices2 = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] >= 1
	end
	c[indices1] .= 2
	c[indices2] .= 3
	c[top7_indices] .= 3
	c .+ 3
end

function axis_change_in_price!(fig, data, impulse; options)
	shocks = impulse_shock(data, impulse)
	model = Model(data, shocks, options)
	sol = solve(model)
	top7_indices = sortperm(shocks.demand_shock, rev = true)[1:7]

	colors = Makie.wong_colors()
	c = get_color(data, shocks)
	ax = Axis(fig[1, 2],
		xlabel = "Change in quantities",
		ylabel = "Change in prices",
		ytickformat = "{:.2f}%",
		xtickformat = "{:.2f}%",
		xlabelsize = 24,
		ylabelsize = 24,
		xticklabelsize = 20,
		yticklabelsize = 20)
	for (i, label) in enumerate(names(impulse)[2:72])
		p = Point2f(100 * (sol.quantities[i] / data.λ[i] - 1), 100 * (sol.prices[i] - 1))
		scatter!(p, markersize = 20, color = colors[c[i]])
		if (i in top7_indices)
			text!(p, text = label, color = :gray70, offset = (0, 10),
				align = (:center, :bottom))
		end
	end
end

function panel(data, impulses; options, name = "panel")
	fig = Figure(size = (1980, 1020))
	axis_change_in_level!(fig, data, impulses, options = options)
	axis_change_in_price!(fig, data, impulses, options = options)
	save("plots/$(name).png", fig)
end

function diff_lambda(data, impulses; options, name = "diff_lambda_imp")
	colors = Makie.wong_colors()
	f = Figure(size = (1980, 1000))
	positions_cge = 4 * range(1, 71) .- 2
	ax = Axis(f[1, 1],
		ytickformat = "{:.2f}%",
		xticks = (positions_cge .+ 1, string.(1:71)),
		xlabelsize = 24,
		ylabelsize = 24,
		xticklabelsize = 20,
		yticklabelsize = 20)
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, options)
	sol = solve(model)
	sol_leontief = solve(Model(data, shocks, Leontief()))

	# Find top 7 shocked sectors
	top7_indices = sortperm(shocks.demand_shock, rev = true)[1:7]

	# Create alpha (transparency) array for highlighting
	alphas = fill(1.0, 71)
	# Highlight top 7 shocked sectors with different alpha
	for i in 1:71
		if !(i in top7_indices)
			alphas[i] = 0.6  # More transparent for non-top sectors
		end
	end

	# Calculate components
	consumption_component = 100 .* data.consumption_share ./ data.λ .* (sol.consumption ./ data.consumption_share .- 1)
	intermediate_component = 100 .* (data.λ .- data.consumption_share) ./ data.λ .* ((sol.quantities - sol.consumption) ./ (data.λ .- data.consumption_share) .- 1)
	leontief_total = 100 .* (sol_leontief.quantities[1:71] ./ data.λ .- 1)

	@info consumption_component[33]
	@info intermediate_component[33]
	positions_leontief = 4 * range(1, 71)
	# Create stacked bars for CGE model
	barplot!(ax,
		positions_cge,
		consumption_component,
		width = 1.8,
		color = [RGBAf(colors[1].r, colors[1].g, colors[1].b, alphas[i]) for i in 1:71],
		label = "CGE Consumption",
	)

	barplot!(ax,
		positions_cge,
		intermediate_component,
		offset = max.(consumption_component, 0.0),
		width = 1.8,
		color = [RGBAf(colors[2].r, colors[2].g, colors[2].b, alphas[i]) for i in 1:71],
		label = "CGE Intermediates",
	)

	# Add Leontief bars
	barplot!(ax,
		positions_leontief,  # Offset to the right
		leontief_total,
		width = 1.8,
		color = [RGBAf(colors[3].r, colors[3].g, colors[3].b, alphas[i]) for i in 1:71],
		label = "Leontief Total",
	)

	# Add sector labels for top 7
	for i in top7_indices
		y_val = leontief_total[i]
		text!(ax, positions_cge[i] + 0.4, y_val + 0.1;
			text = names(impulses)[2:72][i],
			align = (:center, :bottom),
			fontsize = 14,
			rotation = π / 6,  # 30 degrees
		)
	end

	# Create legend
	labels = ["CGE Consumption", "CGE Intermediates", "Leontief Total"]
	elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
	axislegend(ax, elements, labels, position = :rt, labelsize = 38)
	save("plots/$(name).png", f)
end

function effect_of_different_elasticities(shocks, data, gdp_effect_simple; labor_slack_function = full_labor_slack, name = "eg")
	a, b, c, d =
		ThreadsX.map([fill(0.99, 3), fill(0.5, 3), fill(0.2, 3), fill(0.1, 3)]) do elasticity
			BeyondHulten.elasticity_gradient(data, shocks, labor_slack_function, false, elasticity)
		end


    cd_elasticities = CESElasticities(0.99,0.99,0.99)
	sol_cd = solve(Model(data, shocks, CES(cd_elasticities, labor_slack_function)))
	model_leontief = Model(data, shocks, Leontief())
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

function comparison_between_labor_slacks(data, shocks, gdp_effect_simple, title)
	no_ls, ls, ls_emp =
		ThreadsX.map([model -> model.data.labor_share, BeyondHulten.full_labor_slack, BeyondHulten.empircial_labor_slack]) do labor_slack_function
			BeyondHulten.elasticity_gradient(data, shocks, labor_slack_function, false, fill(0.1, 3))
		end

	c = Makie.wong_colors()
	f = Figure(size = (1980, 1000), title = title, color = c)
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.1", xgridvisible = false, titlesize = 30, yticklabelsize = 24, xticklabelsize = 24)

	sol_leontief = solve(Model(data, shocks, Leontief()))

	map_to_gdp(x) = 100 .* reverse(map(x -> real_gdp(x), x))


	for (i, (el, linestyle, slack_type)) in enumerate(zip([no_ls, ls, ls_emp], [:solid, :dash, :dot], ["no labor slack", "calibrated labor slack", "empircial labor slack"]))
		lines!(ax, 0.015 .. 0.9, map_to_gdp(el.ϵ), label = "Elasticity between goods -  $(slack_type)", linewidth = 3, linestyle = linestyle, color = c[1])
		lines!(ax, 0.015 .. 0.9, map_to_gdp(el.θ), label = "Elasticity between labour and goods -  $(slack_type)", linewidth = 3, linestyle = linestyle, color = c[2])
		lines!(ax, 0.015 .. 0.9, map_to_gdp(el.σ), label = "Elasticity of consumption - $(slack_type)", linewidth = 3, linestyle = linestyle, color = c[3])
	end

	lines!(ax, [0.9, 0.015], 100 .* fill(real_gdp(sol_leontief), 2), label = "Leontief model", linewidth = 3, color = c[4])
	lines!(ax, [0.9, 0.015], 100 .* fill(gdp_effect_simple, 2), label = "Initial stimulus", linewidth = 3, color = c[5])
	f[2, 1] = Legend(f, ax, labelsize = 24, tellwidth = false, orientation = :horizontal, nbanks = 4)


	save("plots/comparison_between_labor_slacks.png", f)
end


function labor_slack_gradient(data, impulse)
	shocks = BeyondHulten.impulse_shock(data, impulse)
	gdp_effect_simple = 1 + sum(shocks.demand_shock_raw) ./ sum(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])
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

	f = Figure(size = (1000, 800), color = Makie.wong_colors())
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%", ylabel = "GDP", xlabel = "Labour slack", titlesize = 30, yticklabelsize = 24, xticklabelsize = 24, xlabelsize = 26, ylabelsize = 26)
	lines!(ax, range(100, 0, 100), 100 .* labour_slack_gradient, label = "Real GDP", linewidth = 3)
	lines!(ax, [0, 100], 100 .* fill(real_gdp(sol_leontief), 2), label = "Leontief", linewidth = 3)
	lines!(ax, [0, 100], 100 .* fill(real_gdp(sol_cd), 2), label = "Cobb-Douglas", linewidth = 3)
	lines!(ax, [0, 100], 100 .* fill(gdp_effect_simple, 2), label = "Baseline", linewidth = 3)
	axislegend(ax, position = :rb, labelsize = 28)
	save("plots/labor_slack_gradient.png", f)

	@info real_gdp(sol_cd)
end

function plot_real_gdp_gradient(results; title = "Real GDP", cd, leontief, ylims = (97, 103), initial)
	f = Figure(size = (1980, 1000), title = title, color = Makie.wong_colors())

	ga = f[1, 1] = GridLayout()
	ax = [
		Axis(ga[1, 1], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.9", xgridvisible = false, titlesize = 30, yticklabelsize = 24, xticklabelsize = 24),
		Axis(ga[1, 2], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.5", xgridvisible = false, titlesize = 30, yticklabelsize = 24, xticklabelsize = 24),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.2", xgridvisible = false, titlesize = 30, yticklabelsize = 24, xticklabelsize = 24),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.1", xgridvisible = false, titlesize = 30, yticklabelsize = 24, xticklabelsize = 24)
		]
	ylims!(ax[1], 98, 104)
	ylims!(ax[2], 98, 104)
	ylims!(ax[3], 98, 104)
	ylims!(ax[4], 98, 104)
	map_to_gdp(x) = 100 .* reverse(map(x -> real_gdp(x), x))
	for (i, el) in enumerate(results)
		# Shade area between Leontief and Cobb Douglas
		band!(ax[i], [0.015, 0.9], 100 .* fill(min(leontief, cd), 2), 100 .* fill(max(leontief, cd), 2),
			alpha = 0.2, color = :gray80)

		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.ϵ), label = "Elasticity between goods", linewidth = 3)
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.θ), label = "Elasticity between labour and goods", linewidth = 3)
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.σ), label = "Elasticity of consumption", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(leontief, 2), label = "Leontief model", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(initial, 2), label = "Initial stimulus", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(cd, 2), label = "Cobb Douglas", linewidth = 3)
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 29, tellwidth = false, orientation = :horizontal, nbanks = 2)

	f
end