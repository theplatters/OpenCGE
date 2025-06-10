#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using GLMakie
using StatsBase
using ThreadsX
using ProgressMeter
using Statistics

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
		xtickformat = "{:.2f}%")
	for (i, label) in enumerate(data.io.Sektoren[1:71])
		p = Point2f(100 * (sol.quantities[i] / data.λ[i] - 1), 100 * (sol.prices[i] - 1))
		scatter!(p, markersize = 20, color = colors[c[i]])
		if (i in top7_indices)
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
	ax = Axis(f[1, 1], ytickformat = "{:.2f}%", xticks = ([], []))
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, options)
	sol = solve(model)
	sol_leontief = solve(Model(data, shocks, leontief))

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

	positions_leontief = 4 * range(1, 71)
	positions_cge = 4 * range(1, 71) .- 2
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
		offset = consumption_component,
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
			text = data.io.Sektoren[i],
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

function simulate(shocks, data, gdp_effect_simple; labor_slack_function = full_labor_slack, name = "eg")
	a, b, c, d =
		ThreadsX.map([fill(0.99, 3), fill(0.5, 3), fill(0.2, 3), fill(0.1, 3)]) do elasticity
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

	sol_leontief = solve(Model(data, shocks, leontief))
	multiplier = (1 - sol_leontief.real_gdp) / (1 - gdp_effect_simple)

	@info "Multiplier: ", multiplier
	panel(data, impulses)
	panel(data, impulses, options = ces_options_ls, name = "panel_ls")
	panel(data, impulses, options = ces_options_ls_empirical, name = "panel_ls_empirical")

	simulate(shocks, data, gdp_effect_simple, labor_slack_function = (model -> model.data.labor_share), name = "impulse")
	simulate(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.full_labor_slack, name = "impulse_ls")
	simulate(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.empircial_labor_slack, name = "impulse_ls_empirical")
	diff_lambda(data, impulses)
	diff_lambda(data, impulses, options = ces_options_ls, name = "diff_lambda_imp_ls")
	diff_lambda(data, impulses, options = ces_options_ls_empirical, name = "diff_lambda_imp_ls_empirical")
	labor_slack_gradient(data, impulses)
end


function check_hunch_intermediates_glass(data, impulses)
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, ces_options_ls)
	sol = solve(model)

	f = Figure(size = (1000, 800))
	ax = Axis(f[1, 1], ylabel = "p", xlabel = "Additional labour")

	indices = findall(eachsector(sol)) do x
		index = findfirst(==(x.name), data.io.Sektoren)
		x.price > 1.0 && x.quantity / data.λ[index] < 1
	end

	for (i, label) in enumerate(data.io.Sektoren[1:71])
		p = Point2f(100 * (BeyondHulten.full_labor_slack(model)[i] - data.labor_share[i]), sol.prices[i])
		scatter!(ax, p, markersize = 20)
	end
	f
end
function check_hunch_intermediates_glass_labor(data, impulses)
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, ces_options_ls)
	sol = solve(model)

	sorted_shocks = sortperm(shocks.demand_shock, rev = true)  # Sort indices in descending order
	corr_matrix = cor(data.Ω[sorted_shocks[1:7], 1:71]')

	f = Figure(size = (1200, 800))
	ax = Axis(f[1, 1], xlabel = "Sector", ylabel = "Top 7 Shocked Sectors")

	heatmap!(ax, corr_matrix, colormap = :RdBu, colorrange = (-1, 1))

	# Add sector labels
	ax.xticks = (1:7, data.io.Sektoren[sorted_shocks[1:7]])
	ax.yticks = (1:7, data.io.Sektoren[sorted_shocks[1:7]])
	ax.xticklabelrotation = -π / 2

	Colorbar(f[1, 2], limits = (-1, 1), colormap = :RdBu, label = "Correlation")
	f
	f = Figure(size = (1000, 800))
	ax = Axis(f[1, 1], xlabel = "Sector", ylabel = "Additional labour")
	barplot!(ax, 1:71, BeyondHulten.full_labor_slack(model), color = :steelblue, bar_labels = data.io.Sektoren[1:71],
		label_rotation = -π / 2)

	barplot!(ax, 1:71, data.labor_share, color = :orange)
	f
end

function check_hunch_omega(data, impulses)
	shocks = impulse_shock(data, impulses)
	model = Model(data, shocks, ces_options_ls)
	sol = solve(model)

	f = Figure(size = (1000, 800))
	ax = Axis(f[1, 1], ylabel = "Omega", xlabel = "Sector")
	ax2 = Axis(f[1, 2], ylabel = "Omega", xlabel = "Sector")

	slider = Slider(f[2, 1], range = 1:71, startvalue = 1)
	sector_idx = slider.value

	sector_name = lift(i -> data.io.Sektoren[i], sector_idx)
	omega_values = lift(i -> data.Ω[i, 1:71], sector_idx)
	omega_values2 = lift(i -> data.Ω[1:71, i], sector_idx)

	barplot!(ax, 1:71, omega_values, color = :steelblue)
	barplot!(ax2, 1:71, omega_values2, color = :steelblue)

	ax2.xticks = (1:71, data.io.Sektoren[1:71])
	ax2.xticklabelrotation = -π / 2
	ax2.xticklabelsize = 8
	# Add sector labels on x-axis (rotated for readability)
	ax.xticks = (1:71, data.io.Sektoren[1:71])
	ax.xticklabelrotation = -π / 2
	ax.xticklabelsize = 8

	# Update title with current sector
	ax.ylabel = "Ω value"
	ax.xlabel = "All sectors"
	linkaxes!(ax, ax2)
	# Use on() to update the title when sector_name changes
	on(sector_name) do name
		ax.title = "Selected sector: $name"
	end

	# Set initial title
	ax.title = "Selected sector: $(data.io.Sektoren[1])"

	# Save as HTML
	save("plots/omega_analysis.html", f)
	f
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

function expansion_of_intermediates(data, impulses)
	shocks = impulse_shock(data, impulses)
	model_ls = Model(data, shocks, ces_options_ls)
	sol_ls = solve(model_ls)

	model_ls_emp = Model(data, shocks, ces_options_ls_empirical)
	sol_ls_emp = solve(model_ls_emp)
	(; θ, ϵ, σ) = sol_ls_emp.model.options.elasticities

	q = (data.Ω * sol_ls.prices_raw .^ (1 - θ)) .^ (1 / (1 - θ))
	q_emp = (data.Ω * sol_ls_emp.prices_raw .^ (1 - sol_ls_emp.model.options.elasticities.θ)) .^ (1 / (1 - sol_ls_emp.model.options.elasticities.θ))
	int = ones(71) .^ (-θ) .* (data.Ω' * (ones(71) .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* ones(71) .^ (θ - ϵ) .* (1 .- data.factor_share) .* data.λ))
	int_ls = sol_ls.prices_raw .^ (-θ) .* (data.Ω' * (sol_ls.prices_raw .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* q .^ (θ - ϵ) .* (1 .- data.factor_share) .* sol_ls.quantities))
	int_ls_emp = sol_ls_emp.prices_raw .^ (-θ) .* (data.Ω' * (sol_ls_emp.prices_raw .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* q_emp .^ (θ - ϵ) .* (1 .- data.factor_share) .* sol_ls_emp.quantities))

	@info mean(int_ls ./ int)
	@info mean(int_ls_emp ./ int)

	sol_ls.consumption ./ data.consumption_share + (int_ls ./ (data.λ - data.consumption_share)) - sol_ls.quantities ./ data.λ
end
using Tables, CSV, DataFrames, LinearAlgebra
Omega_R = CSV.File("data/omega_R.csv") |> Tables.matrix
Omega_R = Float64.(Omega_R[[1:71; 73], [2:72; 74]])

consumption_share = data.io[1:length(data.consumption_share), 75] ./ sum(data.io[78, 2:73])
consumption = data.io[1:length(data.consumption_share), 75]

shock = (shocks.demand_shock .- 1) .* consumption

wages = (Vector(data.io[78, 2:72]) ./ data.grossy)

A = vcat(hcat(Matrix(data.io[1:71, 2:72]) ./ (data.grossy'), consumption_share),
	hcat(wages', 0))


q = inv(I - A) * (vcat(shock, 0))
p = ones(length(q))


sol_leontief.quantities[72] ./ sum(mean(col) for col in eachcol(impulses[:, 2:end-2]))
sol_leontief.real_gdp
sum((Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72]) ./ Vector(data.io[findfirst(==("Produktionswert"), data.io.Sektoren), 2:72])) .* (q[1:71])) 
Vector(data.io[findfirst(==("Produktionswert"), data.io.Sektoren), 2:72]) .- data.io[1:71, "Gesamte Verwendung von Gütern"]
sum((Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72]) ./ Vector(data.io[findfirst(==("Produktionswert"), data.io.Sektoren), 2:72])) .* (sol_leontief.quantities[1:71])) 
sum(mean(col) for col in eachcol(impulses[:, 2:end-2]))

writedlm("data/A.csv", A, ',')