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
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Dict{String,Number})
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
  init=[ones(length(model.data.λ)); model.data.λ],
)
  (; data, options, shocks) = model


  #defines the function:
  #defines the concrete problem to be solved (i.e. with inserted parameter values):
  ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
  x = NonlinearSolve.solve(ProbN, reltol=1e-8, abstol=1e-8).u


  n = length(data.factor_share)
  p = x[1:length(data.consumption_share)]
  q = x[(length(data.consumption_share)+1):end]


  labor = options.labor_slack(model)
  (; ϵ, θ, σ) = options.elasticities
  wages = p .* (shocks.supply_shock .^ ((ϵ - 1) / ϵ)) .* (data.factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
  consumption_share = shocks.demand_shock .* data.consumption_share
  inflator = 1 + mean(p, weights(consumption_share)) / mean(wages, weights(data.factor_share))

  consumption = wages' * labor .* consumption_share .* p .^(-σ) 
  df = DataFrames.DataFrame(
    Dict("prices" => p,
      "prices_shifted" => p ./ mean(p, weights(p .* q)),
      "quantities" => q,
      "value_added_relative" => nominal_increase(p, q, model),
      "value_added" => nominal_increase(p, q, model, relative=false),
      "nominal_gdp" => sum(nominal_increase(p, q, model)),
      "real_gdp" => sum(nominal_increase(p, q, model)) / mean(p,weights(data.λ)),
      "real_gdp2a" => sum(nominal_increase(p, q, model)) / mean(wages, weights(q)),
      "real_gdp3" => sum(nominal_increase(p, q, model)) / mean(wages, weights(p .* q)),
      "real_gdp4" => sum(nominal_increase(p, q, model)) / mean(wages, weights(consumption_share)),
      "real_gdp5" => sum(nominal_increase(p, q, model)) / mean(wages, weights(data.factor_share .* p .* q)),
      "real_gdp2" => sum(nominal_increase(p, q, model)) / mean(p, weights(consumption_share)),
      "mean_wages" => mean(p, weights(consumption_share)),
      "real_wage" => mean(wages, weights(data.labor_share)) / mean(p, weights(data.consumption_share)),
      "sectors" => data.io.Sektoren[1:71],
      "consumption" => consumption,
      "wages" => wages),
  )

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
function nominal_gdp(solution::DataFrames.DataFrame; relative=true)
  if !relative
    return solution.nominal_gdp_absolute[1]
  end
  return solution.nominal_gdp[1]
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

function real_gdp(solution::DataFrames.DataFrame; relative=true)
  if !relative
    return solution.real_gdp_absolute[1]
  end
  return solution.real_gdp[1]
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
function nominal_increase(p, q, model::Model; relative=true)
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
function gross_incease(p, q, model; relative=true)
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


