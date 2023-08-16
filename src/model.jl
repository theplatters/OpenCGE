module CESModel
import JuMP
import CSV
import DataFrames

using LinearAlgebra
export solve_ces_model, read_data

struct CESData
  io::DataFrames.DataFrame
  Ω::Matrix{Float64}
  consumption_share::Vector{Float64}
  factor_share::Vector{Float64}
  λ::Vector{Float64}
  labor_share::Vector{Float64}
  A::Vector{Float64}
  B::Vector{Float64}
  θ::Float64
  ϵ::Float64
  σ::Float64
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
  return Ω, consumption_share, factor_share, λ, labor_share
end


function read_data(filename::String)
  filedir = joinpath(pwd(), "data/", filename)
  io = CSV.read(filedir, DataFrames.DataFrame, delim=";", decimal=',', missingstring=["-", "x"])
  DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren)
  io.Sektoren = replace.(io.Sektoren, r"^\s+" => "")
  io = coalesce.(io, 0)

  Ω, consumption_share, factor_share, λ, labor_share = generateData(io)

  return CESData(io, Ω, consumption_share, factor_share, λ, labor_share, ones(71), ones(71), 0.001, 0.5, 0.9)
end

function add_model_variables!(model::JuMP.Model, data::CESData)
  l = length(data.consumption_share)

  JuMP.@variable(model, p[i in 1:l] .>= 0, start = 1)
  JuMP.@variable(model, y[i in 1:l] .>= 0, start = data.λ[i])
end

function add_model_contraints!(model::JuMP.Model, data::CESData)
  p = model[:p]
  y = model[:y]

  JuMP.@NLexpression(model,β[i = 1:71],data.consumption_share[i] * data.B[i] / sum(data.consumption_share[j] * data.B[j] for j in 1:71))
  
  JuMP.@NLexpression(model, q[i=1:71],
    sum(data.Ω[i, j] * p[j]^(1 - data.θ) for j in 1:71)^(1 / (1 - data.θ)))

  JuMP.@NLexpression(model, w[i=1:71],
   p[i] * (data.A[i] ^ ((data.ϵ - 1)/data.ϵ)) * (data.factor_share[i] ^ (1 / data.ϵ)) * (y[i] ^ (1/data.ϵ)) * data.labor_share[i] ^ (-1/data.ϵ))
  JuMP.@NLexpression(model, C, sum(w[j] * data.labor_share[j] for j in 1:71))

  JuMP.@NLconstraint(model, price[i=1:71],
    p[i] == data.A[i]^-1 * (data.factor_share[i] * w[i]^(1 - data.ϵ) + (1 - data.factor_share[i]) * q[i]^(1 - data.ϵ))^(1 / (1 - data.ϵ)))
  
  JuMP.@NLconstraint(model, amount[i=1:71],
    y[i] == p[i]^(-data.θ) * sum(data.A[j]^(data.ϵ - 1) * data.Ω[i, j] * p[j]^data.ϵ * q[j]^(data.θ - data.ϵ) * (1 - data.factor_share[j]) * y[j] for j in 1:71) - p[i]^(-data.σ) * C * β[i])


end



function solve_ces_model(data, optimizer)
  model = JuMP.Model(optimizer)
  data = read_data("I-O_DE2019_formatiert.csv")
  add_model_variables!(model, data)
  add_model_contraints!(model, data)
  p = model[:p]
  JuMP.@objective(model,Min, sum(p[i] - 1 for i in 1:71))
  JuMP.optimize!(model)
  return model
end





end
