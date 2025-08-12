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
	Shocks(supply_shock, demand_shock,zeros(71))
end

function impulse_shock(data, impulses)
	impules_2019_prices = impulses[:, 2:end-2] ./ inflator
	effect = 1 .+  impules_2019_prices ./ data.io[1:71, "Letzte Verwendung von Gütern zusammen"]'
	demand_shock = [mean(col) for col in eachcol(effect[1:2, :])]
	supply_shock = ones(71)
	Shocks(supply_shock, demand_shock, [mean(col) for col in eachcol(impules_2019_prices)])
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


"""
Returns the consumer price index (β ̇ p^(1 - σ))^(1/(1-σ) 

"""
function cpi(sol::Solution)
	σ = sol.model.options.elasticities.σ
	(sol.model.data.consumption_share' * sol.prices_raw .^ (1- σ))^(1/(1 - σ))
end