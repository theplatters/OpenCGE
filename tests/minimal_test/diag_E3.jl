# Test: does an AGGREGATE (uniform) autonomous demand shock engage the
# mobile-labor bridge (eta-dependence), whereas the sectoral shock does not?
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
    labor = sum(L_i) - lb * (w/1.0)^η
    num = cpi - 1.0
    return max(maximum(abs.(zp)), maximum(abs.(mc)), abs(labor), abs(num))
end

# Baseline
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)

# Aggregate (uniform) autonomous demand shock: 0.5% across all sectors.
agg_shocks = Shocks(ones(71), ones(71), zeros(71), fill(0.005, 71), zeros(71))
println("=== AGGREGATE autonomous demand shock (0.5% uniform) ===")
println("eta     maxRes        realGDPΔ%     wage")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    model = mobile_labor_model(data, agg_shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f\n", 0.0, mr, 100*(g/base_gdp - 1), sol.wages[1])
end
# (loop above prints eta=0.0 always due to formatting bug; fix:)
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    model = mobile_labor_model(data, agg_shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol)
    @printf("eta=%6.2f  maxRes=%.2e  realGDPΔ=%+.6f%%  wage=%.6f\n", η, mr, 100*(g/base_gdp - 1), sol.wages[1])
end
println("DONE")
