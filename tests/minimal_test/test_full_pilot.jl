# ═══════════════════════════════════════════════════════════════════════════════
# Full pilot η sweep + variance decomposition — standalone (stdlib only)
# ═══════════════════════════════════════════════════════════════════════════════

using LinearAlgebra
using Printf
using Statistics
using DelimitedFiles

cd("/workspace/BFrep/(3)BeyondHulten")

# ── Data loading (from previous test, proven working) ──
include("data_loader.jl")

data = load_model_data()
(; Ω, consumption_share, factor_share, supply_shock, demand_shock, labor_share, λ, N, grossy, sectors, konsum_col, exporte_col) = data

labor_bar = sum(labor_share)
println("Data loaded: $N sectors, labor_bar = $labor_bar")

# ── Problem function (from previous test, with numeraire) ──
include("problem.jl")

# ── Solver ──
include("solver.jl")

# ═══════════════════════════════════════════════════════════════════════════════
# Task 4: Pilot η sweep
# ═══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("  PILOT η SWEEP")
println("="^70)

η_grid = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
θ_val, ϵ_val, σ_val = 0.5, 0.5, 0.9

params_base = (Ω=Ω, consumption_share=consumption_share, factor_share=factor_share,
               supply_shock=supply_shock, demand_shock=demand_shock, labor_bar=labor_bar,
               θ=θ_val, ϵ=ϵ_val, σ=σ_val)

println("\n  η        Real GDP    Wage      L_total   Conv  Iters")
println("  " * "-"^65)

sweep_sols = []
x0 = [ones(N); λ; 1.0]

for η in η_grid
    params = merge(params_base, (η=η,))
    x_sol, converged, iters = solve_newton(mobile_labor_problem, x0, params, max_iter=500, tol=1e-8)

    p_sol = x_sol[1:N]
    q_sol = x_sol[N+1:2N]
    w_sol = x_sol[2N+1]

    # GDP
    L_i = (p_sol .* (supply_shock .^ ((ϵ_val - 1) / ϵ_val)) .* (factor_share .^ (1 / ϵ_val)) .* (q_sol .^ (1 / ϵ_val)) ./ w_sol) .^ ϵ_val
    L_total = sum(L_i)
    cpi = sum(consumption_share .* p_sol .^ (1 - σ_val))^(1 / (1 - σ_val))
    total_income = w_sol * L_total
    cons_adj = demand_shock .* consumption_share
    cons = total_income .* cons_adj .* (p_sol / cpi) .^ (-σ_val)
    gdp = sum(cons) / sum(consumption_share)

    push!(sweep_sols, (η, gdp, w_sol, L_total, p_sol, q_sol, converged))

    @printf("  %-8.2f %.6f   %.6f   %.6f  %s  %d\n",
        η, gdp, w_sol, L_total, converged ? "✅" : "⚠️", iters)

    global x0 = [p_sol; q_sol; w_sol]
end

# Sectoral analysis
println("\n--- Sectoral Variation Across η ---")
Q_mat = hcat([s[6] for s in sweep_sols]...)
cv_q = vec(std(Q_mat, dims=2))
P_mat = hcat([s[5] for s in sweep_sols]...)

@printf("  Quantity CV: max=%.4f (sector %d: %s), mean=%.4f\n",
    maximum(cv_q), argmax(cv_q), sectors[argmax(cv_q)], mean(cv_q))
@printf("  Price CV:    max=%.4f, mean=%.4f\n", maximum(abs.(std(P_mat, dims=2))), mean(abs.(std(P_mat, dims=2))))
@printf("  Sectors with q-CV > 0.01: %d / %d\n", count(cv_q .> 0.01), N)

# ═══════════════════════════════════════════════════════════════════════════════
# Task 3: Variance Decomposition
# ═══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("  VARIANCE DECOMPOSITION")
println("="^70)

# Full factorial: η × ϵ × θ × σ
η_vals = [0.0, 0.5, 1.0, 5.0, 50.0]     # 5 levels
ϵ_vals = [0.1, 0.5, 0.99]               # 3 levels
θ_vals = [0.1, 0.5, 0.99]               # 3 levels
σ_vals = [0.1, 0.5, 0.99]               # 3 levels
# Total: 5 × 3 × 3 × 3 = 135 evaluations

n_total = length(η_vals) * length(ϵ_vals) * length(θ_vals) * length(σ_vals)
println("\n  Grid: $(length(η_vals))×$(length(ϵ_vals))×$(length(θ_vals))×$(length(σ_vals)) = $n_total evaluations")

# Run the factorial
println("\n  Running (this may take a minute)...")
results = Dict{Tuple{Float64,Float64,Float64,Float64}, Float64}()
eval_count = 0
x0_vd = [ones(N); λ; 1.0]

for η in η_vals, ϵ in ϵ_vals, θ in θ_vals, σ in σ_vals
    global eval_count += 1
    params = (Ω=Ω, consumption_share=consumption_share, factor_share=factor_share,
              supply_shock=supply_shock, demand_shock=demand_shock, labor_bar=labor_bar,
              θ=θ, ϵ=ϵ, σ=σ, η=η)

    try
        x_sol, converged, _ = solve_newton(mobile_labor_problem, x0_vd, params, max_iter=300, tol=1e-7)
        if converged
            p_sol = x_sol[1:N]
            q_sol = x_sol[N+1:2N]
            w_sol = x_sol[2N+1]
            L_i = (p_sol .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (q_sol .^ (1 / ϵ)) ./ w_sol) .^ ϵ
            cpi = sum(consumption_share .* p_sol .^ (1 - σ))^(1 / (1 - σ))
            total_income = w_sol * sum(L_i)
            cons_adj = demand_shock .* consumption_share
            cons = total_income .* cons_adj .* (p_sol / cpi) .^ (-σ)
            gdp = sum(cons) / sum(consumption_share)
            results[(η, ϵ, θ, σ)] = gdp
            global x0_vd = [p_sol; q_sol; w_sol]
        else
            results[(η, ϵ, θ, σ)] = NaN
        end
    catch e
        results[(η, ϵ, θ, σ)] = NaN
    end
end

n_valid = count(!isnan, values(results))
println("  Completed: $n_valid / $n_total converged")

# ── ANOVA decomposition ──
vals = collect(values(results))
valid = .!isnan.(vals)
v = vals[valid]
grand_mean = mean(v)
ss_total = sum((v .- grand_mean) .^ 2)

println("\n  Grand mean GDP: $(round(grand_mean, digits=6))")
println("  SS(total): $(round(ss_total, digits=6))")
println("  GDP range: [$(round(minimum(v), digits=6)), $(round(maximum(v), digits=6))]")

# Compute partial R² for each factor
function partial_r2(factor_name, factor_values, response, grand_mean, ss_total)
    levels = unique(factor_values)
    ss_between = 0.0
    for lvl in levels
        mask = factor_values .== lvl
        group_mean = mean(response[mask])
        n_group = sum(mask)
        ss_between += n_group * (group_mean - grand_mean)^2
    end
    return ss_between / ss_total
end

# Prepare arrays for ANOVA
η_arr = Float64[]
ϵ_arr = Float64[]
θ_arr = Float64[]
σ_arr = Float64[]
gdp_arr = Float64[]

for ((η, ϵ, θ, σ), gdp) in results
    if !isnan(gdp)
        push!(η_arr, η)
        push!(ϵ_arr, ϵ)
        push!(θ_arr, θ)
        push!(σ_arr, σ)
        push!(gdp_arr, gdp)
    end
end

pr2_η = partial_r2("η", η_arr, gdp_arr, grand_mean, ss_total)
pr2_ϵ = partial_r2("ϵ", ϵ_arr, gdp_arr, grand_mean, ss_total)
pr2_θ = partial_r2("θ", θ_arr, gdp_arr, grand_mean, ss_total)
pr2_σ = partial_r2("σ", σ_arr, gdp_arr, grand_mean, ss_total)

println("\n  ┌─────────────────────────────────────┐")
println("  │     PARTIAL R² (Variance Shares)   │")
println("  ├──────────┬───────────┬───────────┤")
println("  │  Factor  │ Partial R²│   % Share │")
println("  ├──────────┼───────────┼───────────┤")
@printf("  │  η       │  %.4f    │  %5.1f %%  │\n", pr2_η, 100*pr2_η/(pr2_η+pr2_ϵ+pr2_θ+pr2_σ))
@printf("  │  ϵ       │  %.4f    │  %5.1f %%  │\n", pr2_ϵ, 100*pr2_ϵ/(pr2_η+pr2_ϵ+pr2_θ+pr2_σ))
@printf("  │  θ       │  %.4f    │  %5.1f %%  │\n", pr2_θ, 100*pr2_θ/(pr2_η+pr2_ϵ+pr2_θ+pr2_σ))
@printf("  │  σ       │  %.4f    │  %5.1f %%  │\n", pr2_σ, 100*pr2_σ/(pr2_η+pr2_ϵ+pr2_θ+pr2_σ))
println("  ├──────────┼───────────┼───────────┤")
total_explained = pr2_η + pr2_ϵ + pr2_θ + pr2_σ
@printf("  │  Total   │  %.4f    │           │\n", total_explained)
@printf("  │  Residual│  %.4f    │           │\n", 1.0 - total_explained)
println("  └──────────┴───────────┴───────────┘")

# ═══════════════════════════════════════════════════════════════════════════════
# Go/No-Go Assessment
# ═══════════════════════════════════════════════════════════════════════════════
println("\n" * "="^70)
println("  GO / NO-GO ASSESSMENT")
println("="^70)

interesting = count(cv_q .> 0.01)
η_dominant = pr2_η > 0.05
max_var = maximum(cv_q)

println("  Sectors with CV > 0.01: $interesting / $N")
println("  Max sectoral CV: $(round(max_var, digits=4))")
println("  η partial R²: $(round(pr2_η, digits=4)) (dominant: $(pr2_η > pr2_ϵ && pr2_η > pr2_θ && pr2_η > pr2_σ))")

go = (interesting >= 5) && (max_var > 0.01)
println("\n  Recommendation: $(go ? "GO ✅" : "NO-GO ⚠️")")
if !go
    if interesting < 5
        println("    - Too few sectors show interesting variation ($interesting < 5)")
    end
    if max_var <= 0.01
        println("    - Maximum sectoral variation too small")
    end
end

println("\n" * "="^70)
println("  DONE")
println("="^70)