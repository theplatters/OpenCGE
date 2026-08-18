# Diag: independently reproduce the second-review's central diagnosis in the
# real 71-sector mobile-labor model:
#   (a) unnormalized demand_shock -> household final demand > income (Walras fail)
#   (b) omitted sector-1 zero-profit residual is materially nonzero
#   (c) goods-market residuals are ~0 (so the omitted eq is the zero-profit one)
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
export variance_decomposition, summary_table, pilot_eta_sweep
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
shocks = standard_shock(data)
model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0)   # η = 0
sol = solve(model)
p, q, w = sol.prices, sol.quantities, sol.wages[1]
(; θ, ϵ, σ, η) = model.options.elasticities
cs = model.data.consumption_share
ds = model.shocks.demand_shock
ss = model.shocks.supply_shock
Ω  = model.data.Ω
fs = model.data.factor_share
labor_bar = model.options.labor_bar

# --- reconstruct final_demand exactly as mobile_labor.jl:150 ---
L_i = sectoral_labor_demand(p, q, w, model)
cpi = sum(cs .* p .^ (1 - σ))^(1 / (1 - σ))
total_income = w * sum(L_i)
final_demand = (total_income .* p .^ (-σ) .* ds .* cs) ./ cpi .^ (-σ)
overspend = sum(final_demand) / total_income

# --- reconstruct sector-1 zero-profit residual (the OMITTED equation) ---
# cost uses SUPPLY_shock (line 156), not demand_shock
intermediate_price = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))
cost = (ss .^ (ϵ - 1) .*
        (fs .* w .^ (1 - ϵ) .+ (1 .- fs) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
zp_residual = p .- cost                     # ~0 for sectors 2..N, large for sector 1
# --- goods-market residuals (all N are enforced; line 153/162) ---
intermediary = p .^ (-θ) .* (Ω' * (p .^ ϵ .* ss .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- fs) .* q))
gm_residual = q .- intermediary .- final_demand
# --- labor-market + numeraire residuals (lines 165/168) ---
labor_res = sum(L_i) - labor_bar * w^η
numeraire_res = cpi - 1.0
total_res_norm = norm([gm_residual; zp_residual; labor_res; numeraire_res])

println("=== Diag: review central diagnosis (71-sector mobile labor, η=0) ===")
println("sum(demand_shock)             = ", round(sum(ds), digits=4),
        "  (level vector: 1.0 baseline + bump in one sector)")
println("total_income                  = ", round(total_income, digits=4))
println("sum(final_demand)             = ", round(sum(final_demand), digits=4))
println("HOUSEHOLD OVERSPEND ratio     = ", round(overspend, digits=4),
        "  (", round(100*(overspend-1), digits=1), "% above income)")
println()
println("max |goods-market residual|   = ", round(maximum(abs.(gm_residual)), digits=6),
        "  (cleared by construction)")
println("zero-profit residual sector1  = ", round(zp_residual[1], digits=4),
        "  (OMITTED eq -> should be ~0 but is not)")
println("max |zero-profit res| sec2..N = ", round(maximum(abs.(zp_residual[2:end])), digits=6))
println("labor-market residual         = ", round(labor_res, digits=6))
println("numeraire (CPI-1) residual    = ", round(numeraire_res, digits=6))
println("TOTAL residual norm           = ", round(total_res_norm, digits=6))
println()
println("=> Household spends ", round(100*(overspend-1), digits=1),
        "% more than income; sector-1 zero-profit residual = ",
        round(zp_residual[1], digits=4), " -> Walras law cannot recover the gap.")
