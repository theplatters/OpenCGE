# Debug: replicate the allocative-wedge math from solved solutions to find
# why the immobile (eta=0) case shows HIGHER real GDP than mobile (eta=1).
using Printf
module GLMakie end
module XLSX end
module BeyondHulten
using NonlinearSolve: NonlinearSolve
using CSV: CSV
using DataFrames
using LineSearches: LineSearches
using LinearAlgebra
using StatsBase
using Printf
using Statistics
using ProgressMeter
using DelimitedFiles
export CESElasticities, MobileLaborCESElasticities, Solution, Shocks, Data, Model
export read_data, standard_shock, solve, CES, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, real_gdp, nominal_gdp, tornqvist_quantity_index
const inflator = 1.46
include("interface.jl")
include("solution.jl")
include("ces.jl")
include("mobile_labor.jl")
include("variance_decomposition.jl")
include("util.jl")
include("impulses.jl")
end
using .BeyondHulten
cd("/workspace/BFrep/(3)BeyondHulten")
data = Data("I-O_DE2019_formatiert.csv")
N = 71
const SEC = 35
const ϵ = 2.0; const θ = 0.5; const σ = 0.9
ss = ones(N); ss[SEC] *= 1.30
shocks = Shocks(ss, ones(N), zeros(N))

for η in [0.0, 1.0]
    model = mobile_labor_model(data, shocks, θ, ϵ, σ, η)
    sol = solve(model)
    p = sol.prices; y = sol.quantities; w = sol.wages[1]
    supply_shock = ss
    factor_share = data.factor_share
    L_opt = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ
    L_base = factor_share .* data.λ
    log_ratio = log.((L_base .+ 1e-12) ./ (L_opt .+ 1e-12))
    mis = 1 - η
    penalty = 0.5 .* factor_share .* (1 .- factor_share) .* ((ϵ - 1) / ϵ) .* (mis .* log_ratio) .^ 2
    alloc_wedge = exp.(-penalty)
    @printf("eta=%.1f  mean(wedge)=%.6f  min(wedge)=%.6f  wedge[35]=%.6f  y[35]=%.6f  sum(y)=%.6f  gdp=%.6f\n",
        η, sum(alloc_wedge)/length(alloc_wedge), minimum(alloc_wedge), alloc_wedge[SEC], y[SEC], sum(y), real_gdp(sol))
end
println("DONE")
