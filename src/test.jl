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
consumption = io[:,"Letzte Verwendung von Gütern zusammen"][1:71]
value_added = Vector(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), 1:72])[2:end]


factor_share = value_added ./ grossy
consumption_share = (I - diagm(1 .- factor_share) * Ω)' * grossy
@views consumption_share[consumption_share.<0] .= 0
consumption_share = consumption_share / sum(consumption_share)
domar_weights = (inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share)
labor = domar_weights .* factor_share

consumption ./ grossy

function problem(X, A, β, Ω, factor_share, ε, θ, σ, labor, B, consumption_share_gross_output)

    N = length(factor_share)
    p = @view X[1:N]
    y = @view X[N+1:end]

    Out = @MVector zeros(eltype(X), 2 * N)

    #Preference Shifter
    β = (B .* β) #/ sum(B .* β)
    #λ = (inv(I - diagm(1 .- α) * Ω)' * β)
    #L = λ .* α
    
    
    labor = inv(I - diagm(1 .- factor_share) * Ω) *  (consumption_share_gross_output .* ((B .* labor) - labor)) + labor

    #labor = B .* y .* α
    #labor = 1.2 * labor
    q = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
    w = p .* (A .^ ((ε - 1) / ε)) .* (factor_share .^ (1 / ε)) .* (y .^ (1 / ε)) .* labor .^ (-1 / ε)
    C = w' * labor

    Out[1:N] = p - (A .^ (ε - 1) .* (factor_share .* w .^ (1 - ε) + (1 .- factor_share) .* q .^ (1 - ε))) .^ (1 / (1 - ε))
    Out[N+1:end] = y - p .^ (-θ) .* (Ω' * (p .^ ε .* A .^ (ε - 1) .* q .^ (θ - ε) .* (1 .- factor_share) .* y)) - C * p .^ (-σ) .* β

    return Out
end

labor = inv(I - diagm(1 .- factor_share) * Ω) *  (consumption_share_gross_output .* ((B .* labor) - labor)) + labor
sum(labor)

consumption_share_gross_output = consumption ./ grossy

f = NonlinearFunction((u, p) -> problem(u, p...))
A = ones(size(grossy))
#A[2] = 1.1

B = ones(71)
names(io)[39]
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

for i in LinRange(1.0, 1.5, 100)

    x0 = x
    B[] = i
    p = [A, consumption_share, Ω, factor_share, ϵ, θ, σ, labor, B]
    ProbN = NonlinearProblem(f, x0, p)

    x = solve(ProbN, NLSolveJL(method=:newton, linesearch=BackTracking()), reltol=1e-8, abstol=1e-8).u

    p = @view x[1:71]
    println("The price after demand dropped to $i is $(p[20])")
    q = @view x[72:end]

    append!(GDP, (p .* (A .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ))' * labor)
end

using Plots
plot(LinRange(1, 1.5, 100), GDP)

#-------------------------------

using JuMP
using Ipopt
include("model.jl")

model = CESModel.solve_ces_model("I-O_DE2019_formatiert.csv",Ipopt.Optimizer)

q = model[:q]
p = model[:p]
y = model[:y]
JuMP.value.(p)