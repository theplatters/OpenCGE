# Milestone D — Honest Variance Decomposition with Sobol indices.
# Reduced grid (5×2×2×2 = 40 pts) to fit memory constraints.
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
export variance_decomposition, summary_table, SobolResult
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

println("═" ^ 72)
println("  Milestone D — Honest Variance Decomposition (Sobol indices)")
println("═" ^ 72)

data = Data("I-O_DE2019_formatiert.csv")

N = 71
SEC = 35
ss = ones(N); ss[SEC] *= 1.30
shocks = Shocks(ss, ones(N), zeros(N))

csv_path = "output/variance_decomposition_sobol.csv"

println("\nReduced grid: η×ε×θ×σ")
println("  η: [0.0, 0.5, 1.0, 2.0]   (4 levels — immobile to mobile)")
println("  ε: [0.5, 2.0]             (2 — complementarity, substitution)")
println("  θ: [0.5, 0.99]            (2 — moderate, high intermed. subst.)")
println("  σ: [0.5, 0.99]            (2 — moderate, high cons. subst.)")
println("  Total: 4×2×2×2 = 32 model evaluations")
println("  Shock: sector-35 (construction) supply +30%")
println()

vd = variance_decomposition(data, shocks;
    η_values = [0.0, 0.5, 1.0, 2.0],
    ϵ_values = [0.5, 2.0],
    θ_values = [0.5, 0.99],
    σ_values = [0.5, 0.99],
    output = :real_gdp, verbose=true)

summary_table(vd; save_csv=csv_path)

println()
println("═" ^ 72)
println("  ASSESSMENT")
println("═" ^ 72)
η_sf = vd.S_f["η"]
η_st = vd.ST_f["η"]
@printf("  η first-order (S_f)   = %.4f  (%.1f%% of variance)\n", η_sf, 100η_sf)
@printf("  η total-order (ST_f)  = %.4f\n", η_st)
@printf("  η interaction share   = %.4f\n", η_st - η_sf)
println()
if vd.n_failed > 0
    println("  ⚠  $(vd.n_failed) grid points failed — Sobol indices may be biased.")
else
    println("  ✅ All $(nrow(vd.grid)) grid points solved successfully.")
end
dom = argmax([vd.S_f[f] for f in vd.factors])
@printf("  Dominant factor: %s (S_f=%.4f)\n", vd.factors[dom], vd.S_f[vd.factors[dom]])
@printf("  Results saved to: %s\n", csv_path)
println("DONE")