# ═══════════════════════════════════════════════════════════════════════════════
# Minimal standalone test — no GLMakie/ModelingToolkit needed
# ═══════════════════════════════════════════════════════════════════════════════
# We include the source files directly, stubbing out the heavy deps.

include("bootstrap.jl")


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
# Test 1: Solve baseline CES model (η=1, cost-minimizing intersectoral allocation)
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
# Test 2: η = 0 (immobile labor — fixed total supply)
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 2: Mobile labor model with η=0.0 (immobile) ---")
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
# Test 3: η = 10 (intersectoral mobility extrapolation)
# ═══════════════════════════════════════════════════════════════════════════════
println("\n--- Test 3: Mobile labor model with η=10.0 (intersectoral mobility extrapolation) ---")
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
