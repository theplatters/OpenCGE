using ProgressMeter, ThreadsX


"""Construct the standard demand shock for one modeled sector."""
function standard_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	demand_shock = ones(n)
	supply_shock = ones(n)
	demand_shock[index] = 1.8097957577943152
	shocks = Shocks(supply_shock, demand_shock, zeros(n))
	return shocks
end


"""
	autonomous_shock(data; sector, autonomous_mult, investment_mult)

Construct additive autonomous and investment-demand multipliers for `sector`.
Inside `MobileLaborCES`, each multiplier is scaled by that sector's consumption
share and baseline aggregate labor income. Other sectors receive zero additive
demand.
"""
function autonomous_shock(data;
	sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten",
	autonomous_mult = 1.8097957577943152,
	investment_mult = 0.0)
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	aut = zeros(n)
	aut[index] = autonomous_mult
	inv = zeros(n)
	inv[index] = investment_mult
	Shocks(ones(n), ones(n), zeros(n), aut, inv)
end


"""Construct the standard technology shock for one modeled sector."""
function standard_tech_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	demand_shock = ones(n)
	supply_shock = ones(n)
	supply_shock[index] = 1.2
	Shocks(supply_shock, demand_shock, zeros(n))
end

"""Construct demand and supply shocks from an impulse-response table."""
function impulse_shock(data, impulses)
	# `impulses` contains year, one column per modeled sector, and wages. Retain every goods
	# sector while excluding the non-sector columns at either end.
	impules_2019_prices = impulses[:, 2:end-1] ./ inflator
	size(impules_2019_prices, 2) == length(data.grossy) ||
		throw(DimensionMismatch("the impulse table must contain one column per modeled sector"))
	n = length(data.grossy)
	last_use = Vector(data.io[1:n, "Letzte Verwendung von Gütern zusammen"])
	effect = 1 .+ impules_2019_prices ./ last_use'
	demand_shock = [mean(col) for col in eachcol(effect[1:min(2, size(effect, 1)), :])]
	supply_shock = ones(n)
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

	return ElasticityGradientSolution(sols_ϵ, sols_θ, sols_σ, labor_reallocation, nominal)
end


"""
Returns the consumer price index (β ̇ p^(1 - σ))^(1/(1-σ) 

"""
function cpi(sol::Solution)
	σ = sol.model.options.elasticities.σ
	(sol.model.data.consumption_share' * sol.prices_raw .^ (1- σ))^(1/(1 - σ))
end
