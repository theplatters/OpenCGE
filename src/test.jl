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
consumption = eachcol(io[:, Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
              sum |>
              x -> getindex(x, 1:71)
value_added = Vector(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), 1:72])[2:end]


α = value_added ./ grossy
consumption_share = (I - diagm(1 .- α) * Ω)' * grossy
@views consumption_share[consumption_share.<0] .= 0
consumption_share = consumption_share / sum(consumption_share)
λ = (inv(I - diagm(1 .- α) * Ω)' * consumption_share)
labor = λ .* α


function problem(X, A, β, Ω, α, ε, θ, σ, L, B)

    N = length(α)
    p = @view X[1:N]
    y = @view X[N+1:end]

    Out = @MVector zeros(eltype(X), 2 * N)

    #Preference Shifter
    β = (B .* β) / sum(B .* β)
    #λ = (inv(I - diagm(1 .- α) * Ω)' * β)
    #L = λ .* α


    q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ε - 1) / ε)) .* (α .^ (1 / ε)) .* (y .^ (1 / ε)) .* L .^ (-1 / ε)
    C = w' * L

    Out[1:N] = p - (A .^ (ε - 1) .* (α .* w .^ (1 - ε) + (1 .- α) .* q .^ (1 - ε))) .^ (1 / (1 - ε))
    Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ε .* A .^ (ε - 1) .* q .^ (θ - ε) .* (1 .- α) .* y)) - C * p .^ (-σ) .* β

    return Out
end

f = NonlinearFunction((u, p) -> problem(u, p...))
A = ones(size(grossy))
A[2] = 1.1

B = ones(71)
B[1:20] .= 1.8
ϵ = 0.5;
θ = 0.001;
σ = 0.9;

x = Complex.([ones(71)..., λ...])
GDP = []

for i in LinRange(1.0, 1.5, 100)

    x0 = x
    B[20] = i
    p = [A, consumption_share, Ω, α, ε, θ, σ, L, B]
    ProbN = NonlinearProblem(f, x0, p)

    x = solve(ProbN, NLSolveJL(method=:newton, linesearch=BackTracking()), reltol=1e-8, abstol=1e-8).u

    p = @view x[1:71]
    println("The price after demand dropped to $i is $(p[20])")
    q = @view x[72:end]

    append!(GDP, (p .* (A .^ ((ε - 1) / ε)) .* (α .^ (1 / ε)) .* (q .^ (1 / ε)) .* L .^ (-1 / ε))' * L)
end


plot(LinRange(1, 0.5, 100), GDP)

#-------------------------------
using JuMP
using Ipopt
include("model.jl")

model = CESModel.solve_ces_model("I-O_DE2019_formatiert.csv",Ipopt.Optimizer)

q = model[:q]
p = model[:p]
y = model[:y]
JuMP.value.(q)