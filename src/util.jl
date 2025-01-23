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


function gradient(data, shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)
	s = copy(sol)
	len = 500
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
				mean_prices[idx] = mean(s.prices_shifted, weights(s.quantities))
			end
		catch e
			@warn e
			gdp[idx] = NaN
			mean_prices[idx] = NaN
			@info idx
		end
	end
	return (gdp, mean_prices)
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
	schedule(t1)
	t2 = @task gradient(data, shocks, labor_slack, labor_reallocation, 2, sol_original, starting_elasticities, nominal)
	schedule(t2)
	t3 = @task gradient(data, shocks, labor_slack, labor_reallocation, 3, sol_original, starting_elasticities, nominal)
	schedule(t3)

	(a, ap) = fetch(t1)
	(b, bp) = fetch(t2)
	(c, cp) = fetch(t3)
	elasticities = starting_elasticities
	return ElasticityGradientSolution(a, b, c, ap, bp, cp, CESElasticities(elasticities...), labor_reallocation, false)
end


function plot_elasticities(results; title = "Real GDP", cd, leontief, ylims = (97, 103))
	f = Figure(size = (1980, 720), title = title)

	ga = f[1, 1] = GridLayout()
	ax = [Axis(ga[1, 1], ylabel = "GDP", ytickformat = "{:.2f}%", title = "0.9"),
		Axis(ga[1, 2], xlabel = "Elasticity", ytickformat = "{:.2f}%", title = "0.5"),
		Axis(ga[2, 1], ytickformat = "{:.2f}%", title = "0.2"),
		Axis(ga[2, 2], ytickformat = "{:.2f}%", title = "0.05")]

	supertitle = Label(f[0, :], title, fontsize = 20, tellwidth = false)
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.σ), label = "Elasticity of consumption")
		lines!(ax[i], [0.9, 0.015], 100 .* fill(leontief, 2), label = "Leontief model", linestyle = :dash)
		#lines!(ax[i], [0.9, 0.015], 100 .* fill(gdp_effect_simple, 2), label = "Baseline Effect", linestyle = :dash)
		lines!(ax[i], [0.9, 0.015], 100 .* fill(cd, 2), label = "Cobb Douglas", linestyle = :dash)
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

	supertitle = Label(f[0, :], title, fontsize = 20, tellwidth = false)
	linkaxes!(ax[1], ax[2], ax[3], ax[4])
	for (i, el) in enumerate(results)
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_ϵ), label = "Elasticity between goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_θ), label = "Elasticity between labour and goods")
		lines!(ax[i], 0.015 .. 0.9, 100 .* reverse(el.mean_prices_σ), label = "Elasticity of consumption")
	end

	f[1, 2] = Legend(f, ax[1], labelsize = 25)

	f
end