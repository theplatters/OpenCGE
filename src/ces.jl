include("interface.jl")


"""
  set_elasticities!(data::CESData, elasticitie::Elasticities)

Set's the elasticities for a given dataset

"""
function set_elasticities!(data, elasticities::AbstractElasticities)
	data.elasticities = elasticities
end


"""
  set_shocks!(data::Data, shocks::Shocks)

Set's the shocks for a given dataset

"""

function set_shocks!(data::AbstractData, shocks::Shocks)
	data.shocks = shocks
end


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
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Vector{<:Number}, sector::Vector{String})

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
function full_demand_labor_allocation(model::Model)
	(; data, shocks, elasticities) = model
	inv(I - diagm(1 .- data.factor_share) * data.Ω) * (data.consumption_share_gross_output .* ((shocks.demand_shock .* data.labor_share) - data.labor_share)) + data.labor_share
end

"""
  problem(X, data::CESData, labor_reallocation)

The objective function as specified in B&F with the added demand shocks, X is the 2*sectors-sized vector,
data contains the parameters and labor_reallocation is a function that specifies how labor is reallocated accross sectors
"""
function problem(X, model::Model{CESElasticities}, labor_reallocation)

	(; data, elasticities, shocks) = model
	N = length(data.factor_share)
	p = @view X[1:N]
	y = @view X[N+1:end]

	(; supply_shock, demand_shock) = shocks
	(; consumption_share, Ω, factor_share) = data
	(; ϵ, θ, σ) = elasticities
	labor = labor_reallocation(model)


	Out = zeros(eltype(X), 2 * N)

	consumption_share = (demand_shock .* consumption_share)

	q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
	w = p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
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
	model::Model{CESElasticities};
	labor_reallocation = full_demand_labor_allocation,
	init = Complex.([ones(71)..., model.data.λ...]),
)
	(; data) = model
	#defines the function:
	f = NonlinearSolve.NonlinearFunction((u, p) -> problem(u, p, labor_reallocation))

	#defines the concrete problem to be solved (i.e. with inserted parameter values):
	ProbN = NonlinearSolve.NonlinearProblem(f, init, model)
	#x_imag gives a vector with prices followed by domar weigths:
	x_imag = NonlinearSolve.solve(ProbN, SciMLNLSolve.NLSolveJL(method = :newton, linesearch = LineSearches.BackTracking()), reltol = 1e-8, abstol = 1e-8).u
	#turns complex numbers into real numbers:
	x = real.(x_imag)
	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]
	df = DataFrames.DataFrame(
		Dict("prices" => p,
			"quantities" => q,
			"value_added_relative" => gross_incease(p, q, model, labor_reallocation),
			"value_added_nominal_relative" => nominal_increase(p, q, model, labor_reallocation),
			"value_added_absolute" => gross_incease(p, q, model, labor_reallocation, relative = false),
			"value_added_nominal_absolute" => nominal_increase(p, q, model, labor_reallocation, relative = false),
			"sectors" => data.io.Sektoren[1:71]),
	)

	return df
end

"""
	nominal_gdp(p,q,data)

Returns the nominal GDP

# Example
```julia-repl
julia> nominal_gdp(ones(76),data.λ,data)
1.0
```
"""
function nominal_gdp(p, q, model)
	(; data, shocks, elasticities) = model
	A = shocks.supply_shock
	(;factor_share, labor_share) = data

	ϵ = elasticities.ϵ
	#put these values into GDP equation from baquee/farhi (difference from real gdp: we do not divide
	# q by prices):
	(p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor_share .^ (-1 / ϵ))' * labor_share
end

"""
	real_gdp(solution::DataFrame;relative)

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
real_gdp(p,q,data, labor_realloc)

Returns the GDP adapted to the pre-shock price level

# Example
```julia-repl
julia> real_gdp(ones(76),data.λ,data,labor_realloc)
1
```
"""
function real_gdp(p, q, model, labor_realloc)
	(; data, shocks, elasticities) = model
	A = shocks.supply_shock
	(;factor_share, labor_share) = data

	ϵ = elasticities.ϵ
	#put these values into real GDP equation from baquee/farhi:
	#q = q ./ p
	((A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor_share .^ (-1 / ϵ))' * labor_share
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
	relative ? sum(solution.value_added_relative) : sum(solution.value_added_absolute)
end
"""
	nominal_increase(q,data)

Returns the increase in output in each sector, in the current price level

# Example
```julia-repl
julia> nominal_increase(data.λ,data,labor_realloc)
0
0
0
```
"""
function nominal_increase(p, q, model, labor_realloc; relative = true)
	(; data, shocks, elasticities) = model
	A = shocks.supply_shock
	(;factor_share, labor_share) = data

	ϵ = elasticities.ϵ
	#put these values into real GDP equation from baquee/farhi:
	#q = q ./ p
	labor = labor_realloc(model)
	w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
	relative ? w .* labor : Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72]) .* w .* labor
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
function gross_incease(p, q, model, labor_realloc; relative = true)
	(; data, shocks, elasticities) = model
	A = shocks.supply_shock
	(;factor_share, labor_share) = data

	ϵ = elasticities.ϵ
	labor = labor_realloc(model)
	w = (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)

	relative ? w .* labor : Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72]) .* w .* labor

end


