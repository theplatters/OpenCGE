# Milestone E verification (Venue 2, then Venue 2+1 layer): an EXPANSIVE
# autonomous export-demand shock should engage the labor channel so that the
# real-GDP response is eta-dependent (eta->infinity > eta=0).
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

# Reconstruct the FULL corrected residual vector (incl. autonomous/investment).
function full_res(model, sol)
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
    debt = sum(p .* total_fd) - total_income   # financing gap = Σ p_i(A+G)
    return (zp=zp, mc=mc, labor=labor, num=num, debt=debt)
end

# True no-shock baseline (all shocks = 1 / 0).
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)

println("=== Milestone E (Venue 2: autonomous export demand) ===")
println("eta     maxRes        realGDPΔ%     wage        debt(financing)")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    shocks = autonomous_shock(data; autonomous_mult = 0.2)   # demand_shock=1, autonomous[constr]=0.2
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_res(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f  %.6f\n",
        η, max(maximum(abs.(r.zp)), maximum(abs.(r.mc)), abs(r.labor), abs(r.num)),
        100*(g/base_gdp - 1), sol.wages[1], r.debt)
end

println("\n=== Milestone E (Venue 2 + Venue 1 layer: debt-financed investment) ===")
for η in [0.0, 1.0, 50.0]
    shocks = autonomous_shock(data; investment_mult = 0.2)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_res(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f  %.6f\n",
        η, max(maximum(abs.(r.zp)), maximum(abs.(r.mc)), abs(r.labor), abs(r.num)),
        100*(g/base_gdp - 1), sol.wages[1], r.debt)
end
println("DONE")
