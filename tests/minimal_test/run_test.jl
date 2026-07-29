# ═══════════════════════════════════════════════════════════════════════════════
# Minimal standalone test — no GLMakie/ModelingToolkit needed
# ═══════════════════════════════════════════════════════════════════════════════
# We include the source files directly, stubbing out the heavy deps.

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

# ── Stub out the modules we don't need (GLMakie, XLSX, etc.) ──
module GLMakie
    # no-op stub
end

module XLSX
    # no-op stub
end

# ── Define the module with only the needed components ──
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

export CESElasticities, LeontiefElasticies, CobbDouglasElasticities, LeontiefElasticiesLabor, MobileLaborCESElasticities
export Solution, eachsector
export Shocks, Data, CBData, read_data_cb
export read_data, calculate_investment!
export gross_increase, nominal_increase, real_gdp, nominal_gdp, full_labor_slack, gdp
export cobb_douglas_costfun, cobb_douglas_consumption, cobb_douglas_intermediary_demand, cobb_douglas_prices!
export solve
export CES, Leontief, CobbDouglas, Solution, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, economy_wide_wage
export elasticities_gradient, standard_shock, standard_tech_shock, ElasticityGradientSolution
export impulse_shock, load_impulses, multiplier, full_labor_slack_alt, inflator, cpi
export plot_real_gdp_gradient, plot_nominal_gdp_gradient, panel, diff_lambda
export comparison_between_labor_slacks, effect_of_different_elasticities
export eta_sweep, eta_sweep_full, EtaSweepResult, real_gdp_sweep, nominal_gdp_sweep
export sectoral_quantities, sectoral_prices
export variance_decomposition, VarianceDecompositionResult, summary_table, pilot_eta_sweep

const inflator = 1.46

include("interface.jl")
include("solution.jl")
include("ces.jl")
include("mobile_labor.jl")
include("variance_decomposition.jl")
include("util.jl")
include("impulses.jl")

end # module

# ═══════════════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════════════
using .BeyondHulten
using Printf

cd("/workspace/BFrep/(3)BeyondHulten")

println("\n" * "="^70)
println("  Testing Mobile Labor CES Model")
println("="^70)

# ── Load data ──
data = Data("I-O_DE2019_formatiert.csv")
N = length(data.factor_share)
println("Data loaded: $N sectors")
println("Sum of labor_share: $(sum(data.labor_share))")
println("Sum of λ: $(sum(data.λ))")

# ── Construct standard shock ──
shocks = standard_shock(data)
println("\nStandard shock constructed")
println("  demand_shock[1:5]: $(shocks.demand_shock[1:5])")
println("  supply_shock[1:5]: $(shocks.supply_shock[1:5])")

# ═══════════════════════════════════════════════════════════════════════════════
# Test 1: Solve baseline CES model (η → ∞, i.e. full labor slack equivalent)
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 1: Mobile labor model with η=1.0 ---")
model1 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)
println("  labor_bar = $(model1.options.labor_bar)")
try
    sol1 = solve(model1)
    println("  ✅ Solved successfully!")
    @printf("  Real GDP:      %.6f\n", real_gdp(sol1))
    @printf("  Nominal GDP:   %.6f\n", nominal_gdp(sol1))
    @printf("  Wage (w):      %.6f\n", sol1.wages[1])
    @printf("  Prices[1:5]:   %s\n", round.(sol1.prices[1:5], digits=4))
    @printf("  Quantities[1:5]: %s\n", round.(sol1.quantities[1:5], digits=4))
    @printf("  Numeraire:     %.6f\n", sol1.numeraire)
catch e
    println("  ❌ Solve failed: $e")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Test 2: η = 0 (perfectly inelastic labor — fixed total supply)
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 2: Mobile labor model with η=0.0 (inelastic) ---")
model2 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0)
try
    sol2 = solve(model2)
    println("  ✅ Solved successfully!")
    @printf("  Real GDP:      %.6f\n", real_gdp(sol2))
    @printf("  Wage (w):      %.6f\n", sol2.wages[1])
catch e
    println("  ❌ Solve failed: $e")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Test 3: η = 10 (highly elastic labor)
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 3: Mobile labor model with η=10.0 (highly elastic) ---")
model3 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 10.0)
try
    sol3 = solve(model3)
    println("  ✅ Solved successfully!")
    @printf("  Real GDP:      %.6f\n", real_gdp(sol3))
    @printf("  Wage (w):      %.6f\n", sol3.wages[1])
catch e
    println("  ❌ Solve failed: $e")
end

# ═══════════════════════════════════════════════════════════════════════════════
# Test 4: η sweep
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 4: η sweep ---")
η_grid = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
try
    sols = eta_sweep(data, shocks, 0.5, 0.5, 0.9, η_grid)
    println("  ✅ Sweep completed!")
    println("\n  η        Real GDP    Δ%")
    println("  " * "-"^40)
    base_gdp = real_gdp(sols[1])
    for (i, η) in enumerate(η_grid)
        @printf("  %-8.2f %.6f  %+.4f%%\n", η, real_gdp(sols[i]), 100*(real_gdp(sols[i]) - base_gdp))
    end

    # Sectoral variation
    Q = sectoral_quantities(EtaSweepResult(η_grid, sols))
    cv = vec(std(Q, dims=2))
    @printf("\n  Sectoral variation (CV across η):\n")
    @printf("    Max:  sector %d = %.4f\n", argmax(cv), maximum(cv))
    @printf("    Mean: %.4f\n", mean(cv))
    @printf("    # sectors with CV > 0.01: %d / %d\n", count(cv .> 0.01), N)
catch e
    println("  ❌ Sweep failed: $e")
    @info "Full error:" e
end

println("\n" * "="^70)
println("  All tests complete")
println("="^70)