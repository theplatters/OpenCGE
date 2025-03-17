function standard_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	shock_amount = 100_000
	demand_shock = ones(71)
	supply_shock = ones(71)
	demand_shock[findfirst(==(sector), data.io.Sektoren)] = 1.8097957577943152
	shocks = Shocks(supply_shock, demand_shock)
	return shocks
end


function standard_tech_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	demand_shock = ones(71)
	supply_shock = ones(71)
	supply_shock[findfirst(==(sector), data.io.Sektoren)] = 1.2
	@info supply_shock[findfirst(==(sector), data.io.Sektoren)]
	Shocks(supply_shock, demand_shock)
end

function impulse_shock(data, impulses)
	effect = 1 .+ impulses[:, 2:end-2] ./ data.io[1:71, "Letzte Verwendung von Gütern zusammen"]'
	demand_shock = [mean(col) for col in eachcol(effect[1:2, :])]
	supply_shock = ones(71)
	Shocks(supply_shock, demand_shock)
end
struct ElasticityGradientSolution
	ϵ::Vector{DataFrame}
	θ::Vector{DataFrame}
	σ::Vector{DataFrame}
	labor_realloc::Bool
	nominal::Bool
end


function gradient(data, shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)
	len = 1000
	sols = Vector{Union{DataFrame, Nothing}}(undef, len)
	arr = copy(el)
	u0 = [ones(71); data.λ]
	@inbounds for (idx, i) in enumerate(range(0.99, 0.015, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		try
			sol_prev = solve(model, init = u0)
			u0 = [sol_prev.prices; sol_prev.quantities]
			sols[idx] = sol_prev
		catch e

			sols[idx] = nothing
			@warn e
			@info idx
		end
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
	t1 = @task gradient(data, shocks, labor_slack, labor_reallocation, 1, sol_original, starting_elasticities, nominal)
	t2 = @task gradient(data, shocks, labor_slack, labor_reallocation, 2, sol_original, starting_elasticities, nominal)
	t3 = @task gradient(data, shocks, labor_slack, labor_reallocation, 3, sol_original, starting_elasticities, nominal)
	schedule(t1)
	schedule(t2)
	schedule(t3)

	sols_ϵ, sols_θ, sols_σ = fetch(t1), fetch(t2), fetch(t3)

	return ElasticityGradientSolution(sols_ϵ, sols_θ, sols_σ, labor_reallocation, nominal)
end


function plot_elasticities(results; title = "Real GDP", cd, leontief, ylims = (97, 103), initial)
	f = Figure(size = (1980, 1000), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "GDP", ytickformat = "{:.2f}%", title = "0.9"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.5"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]

	map_to_gdp(x) = 100 .* reverse(map(x -> x.real_gdp[1], x))
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

function plot_prices(results; title = "Real GDP", cd = sol_cd, ylims = (97, 103))
	f = Figure(size = (1980, 100), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "GDP", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]

	supertitle = Label(f[0, :], title, fontsize = 20, tellwidth = false)
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_σ), label = "Elasticity of consumption")
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


	map_to_real_wages(x) = 100 .* reverse(map(x -> x.real_wage[1], x))
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
	ax = [Axis(ga[1, 1], ylabel = "Wages", ytickformat = "{:.2f}%", title = "0.99"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.7"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]


	map_to_consumption(sols) = 100 .* reverse(map(x -> sum(x.prices .* x.consumption ./ x.numeraire[1]) - x.real_gdp[1], sols))
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, map_to_consumption(el.σ), label = "Elasticity of consumption")
	end

	f[2, 1] = Legend(f, ax[1], labelsize = 25, tellwidth = false, orientation = :horizontal)

	f
end