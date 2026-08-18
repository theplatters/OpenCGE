# Diagnostic: corrected mobile-labor supply anchored at baseline wage w0.
# Reuses the module's data structures and equations, only fixes the
# labor-supply normalization L = Lbar * (w/w0)^eta.

using CSV: CSV
using DataFrames
using LinearAlgebra
using NonlinearSolve: NonlinearSolve
using LineSearches: LineSearches
using StatsBase
using Printf
using Statistics
using ProgressMeter
using DelimitedFiles

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
export sectoral_labor_demand, real_gdp, nominal_gdp
const inflator = 1.46
include("interface.jl")
include("solution.jl")
include("ces.jl")
include("mobile_labor.jl")
include("util.jl")
include("impulses.jl")
end

using .BeyondHulten
using Printf

cd("/workspace/BFrep/(3)BeyondHulten")

data = Data("I-O_DE2019_formatiert.csv")
N = length(data.factor_share)

# Baseline (no shock) shocks to recover baseline wage w0
base_shocks = Shocks(ones(N), ones(N), zeros(N))
m0 = mobile_labor_model(data, base_shocks, 0.5, 0.5, 0.9, 0.0)
sol0 = solve(m0)
w0 = sol0.wages[1]
@printf("Baseline wage w0 = %.6f  (this is what labor supply must be anchored to)\n", w0)

# Corrected problem: labor supply L = Lbar * (w/w0)^eta
function corrected_problem!(out, X, model, w0)
    (; data, options, shocks) = model
    N = length(data.factor_share)
    p = max.(X[1:N], 0.0)
    y = max.(X[N+1:2N], 0.0)
    w = max(X[2N+1], 1e-10)
    (; supply_shock, demand_shock) = shocks
    (; consumption_share, Ω, factor_share) = data
    (; θ, ϵ, σ, η) = options.elasticities

    intermediate_price = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
    L_i = sectoral_labor_demand(p, y, w, model)
    total_income = w * sum(L_i)
    final_demand = (total_income .* (p ./ cpi) .^ (-σ) .* demand_shock .* consumption_share)
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))

    out[1:N-1] .= p[2:N] .- cost[2:N]
    out[N:2N-1] .= y .- intermediary_demand .- final_demand
    out[2N] = sum(L_i) - options.labor_bar * (w / w0)^η
    out[2N+1] = cpi - 1.0
    nothing
end

function corrected_solve(model, w0; init=nothing)
    (; data, options, shocks) = model
    N = length(data.factor_share)
    (; σ) = options.elasticities
    if init === nothing
        init = [ones(N); data.λ; 1.0]
    end
    prob = NonlinearSolve.NonlinearProblem((out,X,_) -> corrected_problem!(out, X, model, w0), init, model)
    x = NonlinearSolve.solve(prob, reltol=1e-8, abstol=1e-8).u
    p = x[1:N]; q = x[N+1:2N]; w = x[2N+1]
    L_i = sectoral_labor_demand(p, q, w, model)
    consumption_share_adj = shocks.demand_shock .* data.consumption_share
    numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
    total_income = w * sum(L_i)
    consumption = total_income .* consumption_share_adj .* (p ./ numeraire) .^ (-σ)
    laspeyres = sum(consumption) / sum(data.consumption_share)
    return (p=p, q=q, w=w, rgdp=laspeyres, L=sum(L_i))
end

shocks = standard_shock(data)
println("\nCorrected η sweep (labor supply anchored at w0):")
println("  η        Real GDP    Total Labor   Wage")
println("  " * "-"^50)
for η in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
    m = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    r = corrected_solve(m, w0)
    @printf("  %-8.2f %.6f    %.6f      %.6f\n", η, r.rgdp, r.L, r.w)
end
