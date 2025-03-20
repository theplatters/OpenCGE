"""
	calculate_investment!(shock::Shocks, data::Data, investment::Number, sector::String)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Vector{<:Number}, sector)

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	for i in eachindex(sector)
		sector_number = findfirst(==(sector[i]), data.io.Sektoren)
		shocks.demand_shock[sector_number] = 1 + investment[i] / consumption[sector_number]
		println("Demand shock to sector $(sector[i]): $(shocks.demand_shock[sector_number])")
	end
end

"""
	calculate_investment!(shock::Shocks, data::Data, investment::Dict)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Dict{String, Number})
	for (sector, investment) in investment
		calculate_investment!(shocks, data, [investment], [sector])
	end
end

"""
	full_labor_slack(model::Model)

Returns the labor vector adjusted, so that labor can be freely reallocated to accomodate for demand shocks
"""
function full_labor_slack(model::Model)
	(; data, shocks) = model
	inv(I - diagm(1 .- data.factor_share) * data.Ω) * (data.consumption_share_gross_output .* ((shocks.demand_shock .* data.labor_share) - data.labor_share)) + data.labor_share
end

function full_labor_slack_constrained(model::Model)
	(; data, shocks) = model
	workforce = 41562 #taken from 
	#https://www-genesis.destatis.de/datenbank/online/table/12211-0001
	total = workforce / (1 - data.unemployment)
	new_labour = inv(I - diagm(1 .- data.factor_share) * data.Ω) * (data.consumption_share_gross_output .* ((shocks.demand_shock .* data.labor_share) - data.labor_share)) + data.labor_share

	if (sum(new_labour) * (1 - data.unemployment) > 1)
		new_labour = inv(1 - data.unemployment) * new_labour / sum(new_labour)
	end
	return new_labour
end

"""
  problem(X, model::Model{CES})

The objective function as specified in B&F with the added demand shocks, X is the 2*sectors-sized vector,
data contains the parameters and labor_reallocation is a function that specifies how labor is reallocated accross sectors
"""
function problem(out::Vector, X::Vector, model::Model{CES})

	(; data, options, shocks) = model
	N = length(data.factor_share)
	p = max.(X[1:N], 0)
	y = max.(X[N+1:end], 0)

	(; supply_shock, demand_shock) = shocks
	(; consumption_share, Ω, factor_share) = data
	(; ϵ, θ, σ) = options.elasticities
	labor = options.labor_slack(model)


	consumption_share = (demand_shock .* consumption_share)

	q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

	w = options.labor_reallocation ?
		ones(Float64, length(p)) :
		p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)

	C = w' * labor

	out[1:N] .= p - (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) + (1 .- factor_share) .* q .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
	out[N+1:end] .= y - p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* q .^ (θ - ϵ) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* consumption_share
	nothing
end


"""
	solve_ces_model(data::CESData, shocks, elasticities,[labor_reallocation, init])
The main function of this module, input the relavant model data, shocks and optionally labor_reallocation and
starting vectors and get back the simulated adapted prices and quantities

"""
function solve(
	model::Model{CES};
	init = [ones(length(model.data.λ)); model.data.λ],
)
	(; data, options, shocks) = model


	#defines the function:
	#defines the concrete problem to be solved (i.e. with inserted parameter values):
	ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
	x = NonlinearSolve.solve(ProbN, reltol = 1e-8, abstol = 1e-8).u


	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]


	labor = options.labor_slack(model)
	(; ϵ, θ, σ) = options.elasticities
	wages = p .* (shocks.supply_shock .^ ((ϵ - 1) / ϵ)) .* (data.factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
	consumption_share = shocks.demand_shock .* data.consumption_share
	consumption = wages' * labor .* consumption_share .* p .^ (-σ)
	laspeyres_index = sum(consumption) / sum(data.consumption_share)
	numeraire = mean(p, weights(data.λ))

	return Solution(p, q, wages, consumption, numeraire, laspeyres_index, (wages' * labor) / numeraire, model)
end

