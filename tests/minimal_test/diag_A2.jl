# Deeper Milestone A check: precision + Tornqvist real-GDP index to see whether
# the mobile-labor bridge produces ANY eta-dependent real response.
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
function full_residuals(model, sol)
    N = length(model.data.factor_share)
    p = sol.prices; q = sol.quantities; w = sol.wages[1]
    (; θ,ϵ,σ,η) = model.options.elasticities
    cs = model.data.consumption_share; Ω=model.data.Ω; fs=model.data.factor_share
    ds = model.shocks.demand_shock; ss = model.shocks.supply_shock; lb = model.options.labor_bar
    intermediate_price = (Ω * p .^ (1-θ)) .^ (1/(1-θ))
    cpi = sum(cs .* p .^ (1-σ))^(1/(1-σ))
    L_i = sectoral_labor_demand(p, q, w, model)
    total_income = w * sum(L_i)
    agg = sum(cs .* ds .* p .^ (1-σ))
    final_demand = (cs .* ds .* total_income .* p .^ (-σ)) ./ agg
    intermediary = p .^ (-θ) .* (Ω' * (p .^ ϵ .* ss .^ (ϵ-1) .* intermediate_price .^ (θ-ϵ) .* (1 .- fs) .* q))
    cost = (ss .^ (ϵ-1) .* (fs .* w .^ (1-ϵ) .+ (1 .- fs) .* intermediate_price .^ (1-ϵ))) .^ (1/(1-ϵ))
    zp = p .- cost; mc = q[2:N] .- intermediary[2:N] .- final_demand[2:N]
    labor = sum(L_i) - lb * (w/1.0)^η; num = cpi - 1.0
    budget = sum(p .* final_demand) - total_income
    return (zp=zp, mc=mc, labor=labor, num=num, budget=budget)
end

noshock = Shocks(ones(71), ones(71), zeros(71))
base_model = mobile_labor_model(data, noshock, 0.5, 0.5, 0.9, 0.0)
base_sol = solve(base_model)
base_q = base_sol.quantities; base_p = base_sol.prices
base_gdp_laspeyres = real_gdp(base_sol)

println("eta     max|zp|      LaspeyresΔ%    Tornqvist(y)Δ%   wage")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    shocks = standard_shock(data)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_residuals(model, sol)
    g_las = 100*(real_gdp(sol)/base_gdp_laspeyres - 1)
    g_tor = 100*(tornqvist_quantity_index(sol.prices, sol.quantities, base_p, base_q) - 1)
    @printf("%6.2f  %.2e  %+.6f%%  %+.6f%%  %.6f\n",
        η, maximum(abs.(r.zp)), g_las, g_tor, sol.wages[1])
end
println("DONE")
