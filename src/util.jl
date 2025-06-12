using ProgressMeter, ThreadsX


function standard_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	shock_amount = 100_000
	demand_shock = ones(71)
	supply_shock = ones(71)
	demand_shock[findfirst(==(sector), data.io.Sektoren)] = 1.8097957577943152
	shocks = Shocks(supply_shock, demand_shock, zeros(71))
	return shocks
end


function standard_tech_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	demand_shock = ones(71)
	supply_shock = ones(71)
	supply_shock[findfirst(==(sector), data.io.Sektoren)] = 1.2
	Shocks(supply_shock, demand_shock)
end

function impulse_shock(data, impulses)
	effect = 1 .+ impulses[:, 2:end-2] ./ data.io[1:71, "Letzte Verwendung von Gütern zusammen"]'
	demand_shock = [mean(col) for col in eachcol(effect[1:2, :])]
	supply_shock = ones(71)
	Shocks(supply_shock, demand_shock, [mean(col) for col in eachcol(impulses[:, 2:end-2])])
end
struct ElasticityGradientSolution
	ϵ::Vector{Solution}
	θ::Vector{Solution}
	σ::Vector{Solution}
	labor_realloc::Bool
	nominal::Bool
end


function gradient(data, shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)::Vector{Solution}
	len = 1000
	sols = Vector{Solution}(undef, len)
	arr = copy(el)
	u0 = [sol.prices; sol.quantities]
	@inbounds for (idx, i) in enumerate(range(0.99, 0.015, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		sol_prev = solve(model, init = u0)
		u0 = [sol_prev.prices_raw; sol_prev.quantities]
		sols[idx] = sol_prev
	end
	return sols
end

function elasticity_gradient(data,
	shocks,
	labor_slack = full_labor_slack,
	labor_reallocation = false,
	starting_elasticities = [0.99, 0.99, 0.99],
	nominal = false,
)


	elasticities = CESElasticities(starting_elasticities...)
	ces = CES(elasticities, labor_slack, labor_reallocation)
	model = Model(data, shocks, ces)
	sol_original = solve(model)
	sols_ϵ, sols_θ, sols_σ = fetch.(Threads.@spawn(gradient($data, $shocks, $labor_slack, $labor_reallocation, i, $sol_original, $starting_elasticities, $nominal)) for i in 1:3)
	writedlm("epsilon_$(starting_elasticities[1]).csv", map(x -> x.real_gdp, sols_ϵ), ',')
	writedlm("theta_$(starting_elasticities[1]).csv", map(x -> x.real_gdp, sols_θ), ',')
	writedlm("sigma_$(starting_elasticities[1]).csv", map(x -> x.real_gdp, sols_σ), ',')

	return ElasticityGradientSolution(sols_ϵ, sols_θ, sols_σ, labor_reallocation, nominal)
end


function plot_real_gdp_gradient(results; title = "Real GDP", cd, leontief, ylims = (97, 103), initial)
	f = Figure(size = (1980, 1000), title = title, color = Makie.wong_colors())

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.9", xgridvisible = false, titlesize = 30),
		Axis(ga[1, 2], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.5", xgridvisible = false, titlesize = 30),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.2", xgridvisible = false, titlesize = 30),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "Developement of GDP with elasticities at 0.1", xgridvisible = false, titlesize = 30)]

	map_to_gdp(x) = 100 .* reverse(map(x -> real_gdp(x), x))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		# Shade area between Leontief and Cobb Douglas
		band!(ax[i], [0.015, 0.9], 100 .* fill(min(leontief, cd), 2), 100 .* fill(max(leontief, cd), 2),
			alpha = 0.2, color = :gray80)

		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.ϵ), label = "Elasticity between goods", linewidth = 3)
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.θ), label = "Elasticity between labour and goods", linewidth = 3)
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.σ), label = "Elasticity of consumption", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(leontief, 2), label = "Leontief model", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(initial, 2), label = "Baseline Effect", linewidth = 3)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(cd, 2), label = "Cobb Douglas", linewidth = 3)
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 29, tellwidth = false, orientation = :horizontal)

	f
end

function plot_nominal_gdp_gradient(results; title = "Nominal GDP", cd, leontief, ylims = (97, 103), initial)
	f = Figure(size = (1980, 1000), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "Nominal GDP", ytickformat = "{:.2f}%", title = "0.9"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.5"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]

	map_to_gdp(x) = 100 .* reverse(map(x -> nominal_gdp(x), x))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_gdp(el.σ), label = "Elasticity of consumption")
		lines!(ax[i], [0.9, 0.015], 100 .* fill(leontief, 2), label = "Leontief model", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(initial, 2), label = "Baseline Effect", linestyle = :dot)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(cd, 2), label = "Cobb Douglas", linestyle = :dash)
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 25, tellwidth = false, orientation = :horizontal)

	f
end

function plot_wages(results; title = "Real Wages", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 1000), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "Wages", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]


	map_to_real_wages(x) = 100 .* reverse(map(x -> mean(x.wages), x))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, map_to_real_wages(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_real_wages(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_real_wages(el.σ), label = "Elasticity of consumption")
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 25, tellwidth = false, orientation = :horizontal)

	f
end

function plot_consumption(results; title = "Real Wages", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 1000), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "Consumption", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]


	map_to_consumption(sols) = 100 .* reverse(map(sol -> sum(sol.consumption), sols))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.σ), label = "Elasticity of consumption")
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 25, tellwidth = false, orientation = :horizontal)

	f
end

function plot_pasche_index(results; title = "Pashe Index", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 1000), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "laspeyres Index", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]


	map_to_lp(sols) = 100 .* reverse(map(x -> x.pashe_index[1], sols))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, map_to_lp(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_lp(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_lp(el.σ), label = "Elasticity of consumption")
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 25, tellwidth = false, orientation = :horizontal)

	f
end
