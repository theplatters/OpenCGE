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
const ces_options_ls_alt = CES(ces_elasticities, BeyondHulten.full_labor_slack_alt, false)
const leontief = Leontief()
const cd_options_ls_empirical = CES(cd_elasticities, BeyondHulten.empircial_labor_slack)
const ces_options_ls_empirical = CES(ces_elasticities, BeyondHulten.empircial_labor_slack)



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

	gdp_effect_simple = 1 + sum(shocks.demand_shock_raw) ./ sum(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])

	panel(data, impulses, options = ces_options)
	panel(data, impulses, options = ces_options_ls, name = "panel_ls")
	panel(data, impulses, options = ces_options_ls_alt, name = "panel_ls_alt")
	panel(data, impulses, options = ces_options_ls_empirical, name = "panel_ls_empirical")

	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = (model -> model.data.labor_share), name = "impulse")
	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.full_labor_slack, name = "impulse_ls")
	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.empircial_labor_slack, name = "impulse_ls_empirical")

	comparison_between_labor_slacks(data, shocks, gdp_effect_simple, "")

	diff_lambda(data, impulses, options = ces_options)
	diff_lambda(data, impulses, options = ces_options_ls, name = "diff_lambda_imp_ls")
	diff_lambda(data, impulses, options = ces_options_ls_alt, name = "diff_lambda_imp_ls_alt")
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


sol_leontief = solve(Model(data, shocks, leontief))
sol = solve(Model(data, shocks, CES(CESElasticities(0.9,0.9,0.9), model -> model.data.labor_share, false)))
sol.real_gdp


consumption_component = 100 .* data.consumption_share ./ data.λ .* (sol.consumption ./ data.consumption_share .- 1)
intermediate_component = 100 .* (data.λ .- data.consumption_share) ./ data.λ .* ((sol.quantities - sol.consumption) ./ (data.λ .- data.consumption_share) .- 1)
leontief_total = 100 .* (sol_leontief.quantities[1:71] ./ data.λ .- 1)


diff = consumption_component + intermediate_component - leontief_total


no_shocks = BeyondHulten.Shocks(ones(71), ones(71), zeros(71))
sol_no_shock = Model(data, no_shocks, ces_options_ls_empirical) |> solve

sol.numeraire
