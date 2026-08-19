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
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)
println("baseline: wage=", base_sol.wages[1], "  sum(L_i)=", sum(sectoral_labor_demand(base_sol.prices, base_sol.quantities, base_sol.wages[1], base_model)), "  labor_bar=", base_model.options.labor_bar)

for m in [0.005, 0.03, 0.10]
    agg = Shocks(ones(71), ones(71), zeros(71), fill(m, 71), zeros(71))
    println("\n=== aggregate autonomous demand shock: $(100m)% uniform ===")
    for η in [1.0, 50.0]
        model = mobile_labor_model(data, agg, 0.5, 0.5, 0.9, η)
        sol = solve(model)
        L = sectoral_labor_demand(sol.prices, sol.quantities, sol.wages[1], model)
        @printf("eta=%6.2f  wage=%.6f  sum(L_i)=%.6f  realGDPΔ=%+.6f%%\n",
            η, sol.wages[1], sum(L), 100*(real_gdp(sol)/base_gdp - 1))
    end
end
println("DONE")
