using CSV
using DataFrames
using MAT, LinearAlgebra, Random, NonlinearSolve
using StaticArrays, Distributed, SciMLNLSolve, LineSearches
using JLD
using Statistics
using SharedArrays

function read_data(filename::String)
    filedir = joinpath(pwd(), "data/", filename)
    io = CSV.read(filedir, DataFrame, delim=";", decimal=',', missingstring=["-", "x"])
    rename!(io, Symbol(names(io)[1]) => :Sektoren)
    io.Sektoren = replace.(io.Sektoren, r"^\s+" => "")
    return io
end

io = read_data("I-O_DE2019_formatiert.csv")


io = coalesce.(io, 0)
Ω = Matrix(io[1:71, 2:72])
Ω = Ω ./ sum(Ω, dims=2)

grossy = io[1:71, "Gesamte Verwendung von Gütern"]
consumption = io[:, "Letzte Verwendung von Gütern zusammen"][1:71]
value_added = Vector(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), 1:72])[2:end]


factor_share = value_added ./ grossy
consumption_share = (I - diagm(1 .- factor_share) * Ω)' * grossy
@views consumption_share[consumption_share.<0] .= 0
consumption_share = consumption_share / sum(consumption_share)
domar_weights = (inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share)
labor = domar_weights .* factor_share
consumption_share_gross_output = consumption ./ grossy

function problem(X, data::CESData)

    N = length(factor_share)
    p = @view X[1:N]
    y = @view X[N+1:end]

    Out = @MVector zeros(eltype(X), 2 * N)




    β = (B .* β)

    labor = min.(1.1 * labor, inv(I - diagm(1 .- factor_share) * Ω) * (consumption_share_gross_output .* ((B .* labor) - labor)) + labor)

    q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ε - 1) / ε)) .* (factor_share .^ (1 / ε)) .* (y .^ (1 / ε)) .* labor .^ (-1 / ε)
    C = w' * labor

    Out[1:N] = p - (A .^ (ε - 1) .* (factor_share .* w .^ (1 - ε) + (1 .- factor_share) .* q .^ (1 - ε))) .^ (1 / (1 - ε))
    Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ε .* A .^ (ε - 1) .* q .^ (θ - ε) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* β

    return Out
end


f = NonlinearFunction((u, p) -> problem(u, p...))
A = ones(size(grossy))

B = ones(71)
B[30:36] .= 1.204
B[39] = 1.204
ϵ = 0.5;
θ = 0.001;
σ = 0.9;
x = Complex.([ones(71)..., domar_weights...])
GDP = []

p = [A, consumption_share, Ω, factor_share, ϵ, θ, σ, labor, B, consumption_share_gross_output]

ProbN = NonlinearProblem(f, x, p)

x = solve(ProbN, NLSolveJL(method=:newton, linesearch=BackTracking()), reltol=1e-8, abstol=1e-8).u

p = @view x[1:71]

q = @view x[72:end]

for i ∈ 1:71
    println("Sector $i $(real.(x[i]))")
end

for i ∈ 72:142
    println("Sector $(i - 71) $(real.(x[i]) - domar_weights[i - 71])")
end


q = q ./ p


(p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ))' * labor


#-------------------------------


include("model.jl")

data = CESModel.read_data("I-O_DE2019_formatiert.csv")

A = ones(71)
A[12] = 0.7
B = ones(71)
B[39] = 1.204
shocks = CESModel.Shocks(A, B)

elasticities = CESModel.Elasticities(0.0001, 0.5, 0.9)

labor_realloc(data) = data.labor_share


p_no_realloc, q_no_realloc = CESModel.solve_ces_model(data, shocks, elasticities, labor_reallocation = labor_realloc)
p,q = CESModel.solve_ces_model(data,shocks,elasticities)

CESModel.real_gdp(p, q, data)
CESModel.nominal_gdp(p, q, data)

CESModel.real_gdp(p_no_realloc, q_no_realloc, data)
CESModel.nominal_gdp(p_no_realloc, q_no_realloc, data)
