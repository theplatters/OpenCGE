module CESModel

import NonlinearSolve
import CSV
import DataFrames
import LineSearches
import SciMLNLSolve
using LinearAlgebra

export Elasticities, Shocks, CESData, solve_ces_model, read_data, calculate_investment!

"""
   Elasticities 

Structure, that stores the 3 elasticities 

θ is the elasticity of production factors
ϵ is the elasticity between work and intermediates
σ the substitution elasticitie

# Examples
```julia-repl
julia> Elasticities(0,5,0.3,0.001) 

```
"""
struct Elasticities
  θ::Float64
  ϵ::Float64
  σ::Float64
end

"""
    Shocks

Two vectors, that have to be the same length as the amount of sectors.
Each entry that differs from one, represents a percentage shock in that sector on demand/supply.

# Example
```julia-repl
julia> Shocks(ones(76),ones(76))
```
"""
struct Shocks
  supply_shock::Vector{Float64}
  demand_shock::Vector{Float64}
end

"""
    CESData
Holds all the relevant data of the problem
"""
mutable struct CESData
  io::DataFrames.DataFrame
  Ω::Matrix{Float64}
  consumption_share::Vector{Float64}
  factor_share::Vector{Float64}
  λ::Vector{Float64}
  labor_share::Vector{Float64}
  consumption_share_gross_output::Vector{Float64}
  shocks::Shocks
  elasticities::Elasticities
  grossy::Vector{Float64}
end
"""
    generateData(io::DataFrames.DataFrame)

Helper function that pulls out the key econometric variables used in the model out of the extended io table
# Example
```julia-repl
julia> Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go = generateData(io)
```
"""
function generateData(io::DataFrames.DataFrame)
  Ω = Matrix(io[1:71, 2:72])
  Ω = Ω ./ sum(Ω, dims=2)

  grossy = io[1:71, "Gesamte Verwendung von Gütern"]
  consumption = eachcol(io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
                sum |>
                x -> getindex(x, 1:71)
  value_added = Vector(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), 1:72])[2:end]


  factor_share = value_added ./ grossy
  consumption_share = (I - diagm(1 .- factor_share) * Ω)' * grossy
  @views consumption_share[consumption_share.<0] .= 0
  consumption_share = consumption_share / sum(consumption_share)
  λ = (inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share)
  labor_share = λ .* factor_share
  consumption_share_gross_output = consumption ./ grossy
  return Ω, consumption_share, factor_share, λ, labor_share, consumption_share_gross_output, grossy
end
"""
  set_elasticities!(data::CESData, elasticitie::Elasticities)

Set's the elasticities for a given dataset

"""
function set_elasticities!(data::CESData, elasticities::Elasticities)
  data.elasticities = elasticities
end


"""
  set_shocks!(data::CESData, shocks::Shocks)

Set's the shocks for a given dataset

"""

function set_shocks!(data::CESData, shocks::Shocks)
  data.shocks = shocks
end


"""
    calculate_investment!(shock::Shocks, data::CESData, investment::Number, sector::Int)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::CESData, investment::Number, sector::Int)
  consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
                sum |>
                x -> getindex(x, 1:71)
  shocks.demand_shock[sector] = 1 + investment / consumption[sector]

end

"""
    calculate_investment!(shock::Shocks, data::CESData, investment::Number, sector::String)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::CESData, investment::Number, sector::String)
  sector_number = findfirst(==(sector), data.io.Sektoren)
  consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
                sum |>
                x -> getindex(x, 1:71)
  shocks.demand_shock[sector_number] = 1 + investment / consumption[sector_number]

end

"""
    read_data(filenem::String)

Given a filename of a IO table located in the /data directory this returns the CESData, where shocks are set to ones
and elasticities are set to the ones presente in the paper by B&F
"""
function read_data(filename::String)
  filedir = joinpath(pwd(), "data/", filename)
  io = CSV.read(filedir, DataFrames.DataFrame, delim=";", decimal=',', missingstring=["-", "x"]) #Read in from CSV
  DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren) #Name the indices after the sectors
  io.Sektoren = replace.(io.Sektoren, r"^\s+" => "") #Remove unneccasary whitespaces
  io = coalesce.(io, 0) #Set NANS to 0

  Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, grossy = generateData(io)

  return CESData(io, Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, Shocks(ones(71), ones(71)), Elasticities(0.0001, 0.5, 0.9), grossy)
end

"""
    full_demand_labor_allocation(data::CESData)

Returns the labor vector adjusted, so that labor can be freely reallocated to accomodate for demand shocks
"""
function full_demand_labor_allocation(data::CESData)
  inv(I - diagm(1 .- data.factor_share) * data.Ω) * (data.consumption_share_gross_output .* ((data.shocks.demand_shock .* data.labor_share) - data.labor_share)) + data.labor_share
end

"""
  problem(X, data::CESData, labor_reallocation)

The objective function as specified in B&F with the added demand shocks, X is the 2*sectors-sized vector,
data contains the parameters and labor_reallocation is a function that specifies how labor is reallocated accross sectors
"""
function problem(X, data::CESData, labor_reallocation)

  N = length(data.factor_share)
  p = @view X[1:N]
  y = @view X[N+1:end]

  A = data.shocks.supply_shock
  B = data.shocks.demand_shock
  β = data.consumption_share
  Ω = data.Ω
  factor_share = data.factor_share
  labor = labor_reallocation(data)

  ϵ = data.elasticities.ϵ
  θ = data.elasticities.θ
  σ = data.elasticities.σ

  Out = zeros(eltype(X), 2 * N)

  β = (B .* β)

  q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
  w = p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
  C = w' * labor

  Out[1:N] = p - (A .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) + (1 .- factor_share) .* q .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
  Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ϵ .* A .^ (ϵ - 1) .* q .^ (θ - ϵ) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* β

  return Out
end
"""
    solve_ces_model(data::CESData, shocks, elasticities,[labor_reallocation, init])
The main function of this module, input the relavant model data, shocks and optionally labor_reallocation and
starting vectors and get back the simulated adapted prices and quantities

"""
function solve_ces_model(
  data::CESData,
  shocks, elasticities;
  labor_reallocation=full_demand_labor_allocation,
  init=Complex.([ones(71)..., data.λ...])
)
  set_elasticities!(data, elasticities)
  set_shocks!(data, shocks)
  f = NonlinearSolve.NonlinearFunction((u, p) -> problem(u, p, labor_reallocation))



  ProbN = NonlinearSolve.NonlinearProblem(f, init, data)

  x_imag = NonlinearSolve.solve(ProbN, SciMLNLSolve.NLSolveJL(method=:newton, linesearch=LineSearches.BackTracking()), reltol=1e-8, abstol=1e-8).u
  x = real.(x_imag)
  p = @view x[1:length(data.consumption_share)]
  q = @view x[(length(data.consumption_share)+1):end]
  return p, q
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
function nominal_gdp(p, q, data)
  A = data.shocks.supply_shock
  factor_share = data.factor_share
  labor_share = data.labor_share

  ϵ = data.elasticities.ϵ

  (p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor_share .^ (-1 / ϵ))' * labor_share
end

"""
    real_gdp(p,q,data)

Returns the GDP adapted to the pre-shock price level

# Example
```julia-repl
julia> real_gdp(ones(76),data.λ,data)
1
```
"""
function real_gdp(p, q, data)
  A = data.shocks.supply_shock
  factor_share = data.factor_share
  labor_share = data.labor_share
  ϵ = data.elasticities.ϵ

  q = q ./ p
  (p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor_share .^ (-1 / ϵ))' * labor_share
end

"""
    nominal_increase(q,data)

Returns the increase in output in each sector, in the current price level

# Example
```julia-repl
julia> nominal_increase(data.λ,data)
0
0
0
```
"""
function nominal_increase(q, data)
  return (q - data.λ) * data.grossy
end

"""
    gross_increase(p,q,data)

Returns the increase in output in each sector, in the price level before the shock

# Example
```julia-repl
julia> gross_incease(ones(76),data.λ,data)
zeros(76)
```
"""
function gross_incease(p, q, data)
  return (q ./ p - data.λ) * data.grossy
end

end