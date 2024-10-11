include("interface.jl")




"""
	calculate_investment!(shock::Shocks, data::CESData, investment::Number, sector::Int)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Vector{<:Number}, sector::Vector{Int})

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	for i in 1:length(sector)
		shocks.demand_shock[sector[i]] = 1 + investment[i] / consumption[sector[i]]
		println("Demand shock to sector $(sector): $(shocks.demand_shock[sector[i]])")
	end

end

"""
	calculate_investment!(shock::Shocks, data::Data, investment::Number, sector::String)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Vector{<:Number}, sector)

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	for i in 1:length(sector)
		sector_number = findfirst(==(sector[i]), data.io.Sektoren)
		shocks.demand_shock[sector_number] = 1 + investment[i] / consumption[sector_number]
		println("Demand shock to sector $(sector[i]): $(shocks.demand_shock[sector_number])")
	end

end


"""
	full_demand_labor_allocation(data::CESData)

Returns the labor vector adjusted, so that labor can be freely reallocated to accomodate for demand shocks
"""
function full_labor_slack(model::Model)
	(; data, shocks, options) = model
	elasticites = options.elasticities
	inv(I - diagm(1 .- data.factor_share) * data.Ω) * (data.consumption_share_gross_output .* ((shocks.demand_shock .* data.labor_share) - data.labor_share)) + data.labor_share
end

"""
  problem(X, model::Model{CES})

The objective function as specified in B&F with the added demand shocks, X is the 2*sectors-sized vector,
data contains the parameters and labor_reallocation is a function that specifies how labor is reallocated accross sectors
"""
function problem(X, model::Model{CES})

	(; data, options, shocks) = model
	N = length(data.factor_share)
	p = @view X[1:N]
	y = @view X[N+1:end]

	(; supply_shock, demand_shock) = shocks
	(; consumption_share, Ω, factor_share) = data
	(; ϵ, θ, σ) = options.elasticities
	labor = options.labor_slack(model)


	Out = zeros(eltype(X), 2 * N)

	consumption_share = (demand_shock .* consumption_share)

	q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

	w = options.labor_reallocation ?
		ones(Float64, length(p)) :
		p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)

	C = w' * labor

	Out[1:N] = p - (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) + (1 .- factor_share) .* q .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
	Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* q .^ (θ - ϵ) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* consumption_share

	return Out
end

"""
	solve_ces_model(data::CESData, shocks, elasticities,[labor_reallocation, init])
The main function of this module, input the relavant model data, shocks and optionally labor_reallocation and
starting vectors and get back the simulated adapted prices and quantities

"""
function solve(
	model::Model{CES};
	init = Complex.([ones(71)..., model.data.λ...]),
)
	(; data, options,shocks) = model
	#defines the function:
	f = NonlinearSolve.NonlinearFunction((u, p) -> problem(u, p))

	#defines the concrete problem to be solved (i.e. with inserted parameter values):
	ProbN = NonlinearSolve.NonlinearProblem(f, init, model)
	#x_imag gives a vector with prices followed by domar weigths:
	x_imag = NonlinearSolve.solve(ProbN, SciMLNLSolve.NLSolveJL(method = :newton, linesearch = LineSearches.BackTracking()), reltol = 1e-8, abstol = 1e-8).u
	#turns complex numbers into real numbers:
	x = real.(x_imag)
	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]

	if (options.labor_reallocation)
		df = DataFrames.DataFrame(
			Dict("prices" => p,
				"quantities" => q,
				"sectors" => data.io.Sektoren[1:71],
				"gdp" => ((shocks.demand_shock .* data.consumption_share)' * p .^ (1 - options.elasticities.σ)) ^ (1 / (options.elasticities.σ - 1)),
			))

	else
		df = DataFrames.DataFrame(
			Dict("prices" => p,
				"quantities" => q,
				"value_added_relative" => gross_incease(p, q, model),
				"value_added_nominal_relative" => nominal_increase(p, q, model),
				"value_added_absolute" => gross_incease(p, q, model, relative = false),
				"value_added_nominal_absolute" => nominal_increase(p, q, model, relative = false),
				"sectors" => data.io.Sektoren[1:71]),
		)
	end

	return df
end


"""
	nominal_gdp(solution::DataFrame; relative)

Returns the nominal GDP relative to the preshockgdp, or in € if relative = true is set

```julia-repl
julia> nominal_gdp(sol,relative = true)
  112412.01
julia> nominal_gdp(sol)
  1.023
```
"""
function nominal_gdp(solution::DataFrames.DataFrame; relative = true)
	relative ? sum(solution.value_added_nominal_relative) : sum(solution.value_added_nominal_absolute)
end


"""
	real_gdp(solution::DataFrame;relative)

Returns the GDP relative to the preshockgdp, or in € if relative = true is set

```julia-repl
julia> real_gdp(sol,relative = true)
  112412.01
julia> real_gdp(sol)
```
"""

function real_gdp(solution::DataFrames.DataFrame; relative = true)
	if hasproperty(solution, :gdp)
		return solution.gdp[1]
	end
	relative ? sum(solution.value_added_relative) : sum(solution.value_added_absolute)
end
"""
	nominal_increase(p, q, model)

Returns the increase in output in each sector, in the current price level

# Example
```julia-repl
julia> nominal_increase(data.λ,data,labor_realloc)
0
0
0
```
"""
function nominal_increase(p, q, model::Model; relative = true)
	(; data, shocks, options) = model
	A = shocks.supply_shock
	(; factor_share) = data

	(; ϵ, σ) = options.elasticities
	#put these values into real GDP equation from baquee/farhi:
	#q = q ./ p
	labor = options.labor_slack(model)

	w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
	relative ? w .* labor : sum(Vector(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])) * w .* labor
end


"""
	gross_increase(p,q,data, labor_realloc; relative)

Returns the increase in output in each sector, in the price level before the shock. If relative

# Example
```julia-repl
julia> gross_incease(ones(76),data.λ,data)
zeros(76)
```
"""
function gross_incease(p, q, model; relative = true)
	(; data, shocks, options) = model
	A = shocks.supply_shock
	(; factor_share) = data

	(; ϵ, σ) = options.elasticities
	labor = options.labor_slack(model)
	w = options.labor_reallocation ?
		ones(length(p)) :
		(A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)


	relative ? w .* labor : sum(Vector(data.io[1:71, "Letzte Verwendung von Gütern zusammen"])) * w .* labor

end


