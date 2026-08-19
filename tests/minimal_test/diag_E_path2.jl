# Path 2 verification: sticky-wage / unemployment closure.
# Compare :mobile (wage market-clears, wage pinned -> employment fixed) vs
# :fixed (wage sticky at 1 -> employment absorbs the shock).
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
export read_data, standard_shock, autonomous_shock, solve, CES, MobileLaborCES, mobile_labor_model
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

function maxres(model, sol)
    N = length(model.data.factor_share)
    p = sol.prices; q = sol.quantities; w = sol.wages[1]
    (; θ,ϵ,σ,η) = model.options.elasticities
    cs = model.data.consumption_share; Ω=model.data.Ω; fs=model.data.factor_share
    ds = model.shocks.demand_shock; ss = model.shocks.supply_shock; lb = model.options.labor_bar
    ax = model.shocks.autonomous_demand; gv = model.shocks.investment_shock
    intermediate_price = (Ω * p .^ (1-θ)) .^ (1/(1-θ))
    cpi = sum(cs .* p .^ (1-σ))^(1/(1-σ))
    L_i = sectoral_labor_demand(p, q, w, model)
    total_income = w * sum(L_i)
    agg = sum(cs .* ds .* p .^ (1-σ))
    c = (cs .* ds .* total_income .* p .^ (-σ)) ./ agg
    cons_base = sum(data.labor_share)
    A = ax .* data.consumption_share .* cons_base; G = gv .* data.consumption_share .* cons_base
    total_fd = c .+ A .+ G
    intermediary = p .^ (-θ) .* (Ω' * (p .^ ϵ .* ss .^ (ϵ-1) .* intermediate_price .^ (θ-ϵ) .* (1 .- fs) .* q))
    cost = (ss .^ (ϵ-1) .* (fs .* w .^ (1-ϵ) .+ (1 .- fs) .* intermediate_price .^ (1-ϵ))) .^ (1/(1-ϵ))
    zp = p .- cost
    mc = q[1:N-1] .- intermediary[1:N-1] .- total_fd[1:N-1]
    if model.options.closure == :fixed
        labor = w - 1.0
    else
        labor = sum(L_i) - lb * (w/1.0)^η
    end
    num = cpi - 1.0
    return max(maximum(abs.(zp)), maximum(abs.(mc)), abs(labor), abs(num))
end

# Baseline (no shock)
base = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.5), data))
base_sol = solve(base)
base_gdp = real_gdp(base_sol)
Lbar = base.options.labor_bar

function run(closure, mult)
    shocks = autonomous_shock(data; autonomous_mult = mult)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.5; closure = closure)
    sol = solve(model)
    mr = maxres(model, sol)
    L = sum(sectoral_labor_demand(sol.prices, sol.quantities, sol.wages[1], model))
    @printf("%-7s  maxRes=%.2e  realGDPΔ=%+.6f%%  emp=%.6f  unemp=%.6f  wage=%.6f\n",
        closure, mr, 100*(real_gdp(sol)/base_gdp - 1), L, Lbar - L, sol.wages[1])
end

println("Construction autonomous shock (mult=0.2). Baseline emp=", Lbar, " wage=1.0")
print("  :mobile  "); run(:mobile, 0.2)
print("  :fixed   "); run(:fixed, 0.2)
println("\nConstruction autonomous shock (mult=0.5)")
print("  :mobile  "); run(:mobile, 0.5)
print("  :fixed   "); run(:fixed, 0.5)
println("DONE")
