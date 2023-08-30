module CESModel

import NonlinearSolve
import CSV
import DataFrames
import LineSearches
import SciMLNLSolve 
using LinearAlgebra

export Elasticities, Shocks, CESData, solve_ces_model, read_data


struct Elasticities
  θ::Float64
  ϵ::Float64
  σ::Float64
end

struct Shocks
  supply_shock::Vector{Float64}
  demand_shock::Vector{Float64}
end

struct CESData
  io::DataFrames.DataFrame
  Ω::Matrix{Float64}
  consumption_share::Vector{Float64}
  factor_share::Vector{Float64}
  λ::Vector{Float64}
  labor_share::Vector{Float64}
  consumption_share_gross_output::Vector{Float64}
  shocks::Shocks
  elasticities::Elasticities
end

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
  return Ω, consumption_share, factor_share, λ, labor_share, consumption_share_gross_output
end

function set_elasticities!(data::CESData, elasticities::Elasticities)
  data.elasticities = elasticities
end

function set_shocks!(data::CESData, shocks::Shocks)
  data.shocks = shocks
end

function read_data(filename::String)
  filedir = joinpath(pwd(), "data/", filename)
  io = CSV.read(filedir, DataFrames.DataFrame, delim=";", decimal=',', missingstring=["-", "x"])
  DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren)
  io.Sektoren = replace.(io.Sektoren, r"^\s+" => "")
  io = coalesce.(io, 0)

  Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go = generateData(io)

  return CESData(io, Ω, consumption_share, factor_share, λ, labor_share,consumption_share_go, Shocks(ones(71),ones(71)), Elasticities(0.0001,0.5,0.9))
end



function problem(X, data::CESData)

  N = length(factor_share)
  p = @view X[1:N]
  y = @view X[N+1:end]

  A = data.shoks.supply_shock
  B = data.shoks.demand_shock
  β = data.consumption_share
  consumption_share_gross_output = data.consumption_share_gross_output
  Ω = data.Ω
  factor_share = data.factor_share
  labor = data.labor

  ϵ = data.elasticities.ϵ
  θ = data.elasticities.θ
  σ = data.elasticities.σ

  Out = zeros(eltype(X), 2 * N)

  β = (B .* β)


  labor = min.(1.1 * labor,
    inv(I - diagm(1 .- factor_share) * Ω) * (consumption_share_gross_output .* ((B .* labor) - labor)) + labor)

  q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
  w = p .* (A .^ ((ε - 1) / ε)) .* (factor_share .^ (1 / ε)) .* (y .^ (1 / ε)) .* labor .^ (-1 / ε)
  C = w' * labor

  Out[1:N] = p - (A .^ (ε - 1) .* (factor_share .* w .^ (1 - ε) + (1 .- factor_share) .* q .^ (1 - ε))) .^ (1 / (1 - ε))
  Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ε .* A .^ (ε - 1) .* q .^ (θ - ε) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* β

  return Out
end

function solve_ces_model(data, shocks, elasticities)
  set_elasticities!(data, elasticities)
  set_shocks!(data, shocks)
  f = NonlinearSolve.NonlinearFunction((u, p) -> problem(u, p...))

  ProbN = NonlinearSolve.NonlinearProblem(f, x, data)

  x = NonlinearSolve.solve(ProbN, SciMLNLSolve.NLSolveJL(method=:newton, linesearch=LineSearches.BackTracking()), reltol=1e-8, abstol=1e-8).u

  p = @view x[1:71]
  q = @view x[72:end]
  return p, q
end

function nominal_gdp(p,q,data)
  A = data.shoks.supply_shock
  factor_share = data.factor_share
  labor = data.labor


  (p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ))' * labor
end

function real_gdp(p,q,data)
  A = data.shoks.supply_shock
  factor_share = data.factor_share
  labor = data.labor

  q = q ./ p
  (p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ))' * labor
end




end
