# ═══════════════════════════════════════════════════════════════════════════════
# Variance Decomposition & η Sweep Infrastructure
# ═══════════════════════════════════════════════════════════════════════════════
#
# Implements:
#   1. η sweep — run the mobile-labor CES model across a grid of η values
#   2. Variance decomposition — factorial design over (η, ε, θ, σ) with
#      partial R² to quantify each elasticity's contribution to output variance
#
# The variance decomposition answers: "How much of the variation in GDP (or
# sectoral output) is driven by the labor supply elasticity η vs. the
# technology elasticities (ε, θ, σ)?"
#
# Method: ANOVA-style decomposition on a full factorial grid.
#   - For each factor f, partial R²(f) = SS(f) / SS(total)
#   - SS(f) = sum over all levels of f of [mean(y|f=level) - mean(y)]² · n_level
#   - This gives the share of total variance attributable to each factor.
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

using ProgressMeter
using Printf
using Statistics

"""
    eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=nothing)

Run the mobile-labor CES model across a vector of η values.
Returns a vector of `Solution` objects.

# Example
```julia
η_grid = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, Inf]
sols = eta_sweep(data, shocks, 0.5, 0.5, 0.9, η_grid)
```
"""
function eta_sweep(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η_values::Vector{Float64}; labor_bar::Union{Float64, Nothing}=nothing)
    N = length(data.factor_share)
    results = Vector{Solution}(undef, length(η_values))

    # Use previous solution as warm start
    init = [ones(N); data.λ; 1.0]

    @showprogress "η sweep: " for (idx, η) in enumerate(η_values)
        model = mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=labor_bar)
        sol = solve(model; init=init)
        results[idx] = sol
        # Warm start next iteration
        init = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]
    end

    return results
end

"""
    EtaSweepResult

Container for η sweep results with convenience accessors.
"""
struct EtaSweepResult
    η_values::Vector{Float64}
    solutions::Vector{Solution}
end

function Base.getindex(esr::EtaSweepResult, i::Int)
    return esr.η_values[i], esr.solutions[i]
end

"""
    real_gdp_sweep(esr::EtaSweepResult)

Extract real GDP values from an η sweep.
"""
function real_gdp_sweep(esr::EtaSweepResult)
    [real_gdp(sol) for sol in esr.solutions]
end

"""
    nominal_gdp_sweep(esr::EtaSweepResult)
"""
function nominal_gdp_sweep(esr::EtaSweepResult)
    [nominal_gdp(sol) for sol in esr.solutions]
end

"""
    sectoral_quantities(esr::EtaSweepResult)

Extract a matrix of sectoral quantities (sectors × η values).
"""
function sectoral_quantities(esr::EtaSweepResult)
    N = length(esr.solutions[1].quantities)
    M = length(esr.solutions)
    Q = zeros(N, M)
    for (j, sol) in enumerate(esr.solutions)
        Q[:, j] = sol.quantities
    end
    return Q
end

"""
    sectoral_prices(esr::EtaSweepResult)

Extract a matrix of sectoral prices (sectors × η values).
"""
function sectoral_prices(esr::EtaSweepResult)
    N = length(esr.solutions[1].prices)
    M = length(esr.solutions)
    P = zeros(N, M)
    for (j, sol) in enumerate(esr.solutions)
        P[:, j] = sol.prices
    end
    return P
end

# ═══════════════════════════════════════════════════════════════════════════════
# Variance Decomposition
# ═══════════════════════════════════════════════════════════════════════════════

"""
    VarianceDecompositionResult

Results of a variance decomposition.

- `factors`: names of the factors (e.g., ["η", "ϵ", "θ", "σ"])
- `partial_r2`: partial R² for each factor (sums to ≤ 1; residual = interactions)
- `grid`: the factorial grid used
- `values`: the output values at each grid point
- `output_name`: name of the output variable decomposed
"""
struct VarianceDecompositionResult
    factors::Vector{String}
    partial_r2::Dict{String, Float64}
    grid::DataFrame
    values::Vector{Float64}
    output_name::String
end

"""
    variance_decomposition(data, shocks; η_values, ϵ_values, θ_values, σ_values,
                           output=:real_gdp, labor_bar=nothing)

Run a full factorial design over (η, ε, θ, σ) and decompose the variance
of the output variable into contributions from each elasticity.

# Arguments
- `data::Data`: calibration data
- `shocks::Shocks`: demand/supply shocks
- `η_values`, `ϵ_values`, `θ_values`, `σ_values`: grids for each elasticity
- `output`: Symbol selecting the output variable. Options:
    - `:real_gdp` — real GDP index
    - `:nominal_gdp` — nominal GDP
    - `:sectoral_q` — sectoral quantities (returns per-sector decomposition)
    - `:sectoral_p` — sectoral prices (returns per-sector decomposition)

# Returns
- `VarianceDecompositionResult` with partial R² for each factor

# Method
Uses ANOVA-style one-at-a-time decomposition:
    SS(f) = Σ_levels [mean(y | f=level) - grand_mean]² × n_level
    partial_R²(f) = SS(f) / SS_total

This captures the main effects. Interactions are in the residual.
"""
function variance_decomposition(
    data::Data,
    shocks::Shocks;
    η_values::Vector{Float64} = [0.0, 0.5, 1.0, 2.0, 10.0],
    ϵ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    θ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    σ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    output::Symbol = :real_gdp,
    labor_bar::Union{Float64, Nothing} = nothing,
    verbose::Bool = true,
)

    # Build the full factorial grid
    grid = DataFrame()
    grid.η = Float64[]
    grid.ϵ = Float64[]
    grid.θ = Float64[]
    grid.σ = Float64[]

    for η in η_values, ϵ in ϵ_values, θ in θ_values, σ in σ_values
        push!(grid, (η, ϵ, θ, σ))
    end

    n = nrow(grid)
    if verbose
        @info "Variance decomposition: $n model evaluations across $(length(η_values))×$(length(ϵ_values))×$(length(θ_values))×$(length(σ_values)) grid"
    end

    # ── Run the model for each grid point ──
    N = length(data.factor_share)
    values = Vector{Float64}(undef, n)
    sectoral_values = output in (:sectoral_q, :sectoral_p) ? zeros(N, n) : nothing

    init = [ones(N); data.λ; 1.0]

    @showprogress "Variance decomposition: " for i in 1:n
        η = grid.η[i]
        ϵ = grid.ϵ[i]
        θ = grid.θ[i]
        σ = grid.σ[i]

        # Handle η = Inf (perfectly elastic labor → use full labor slack logic)
        if isinf(η)
            # For η → ∞, w is pinned at the reservation wage (w=1 in real terms)
            # We approximate with a large finite η
            η_eff = 1e6
        else
            η_eff = η
        end

        model = mobile_labor_model(data, shocks, θ, ϵ, σ, η_eff; labor_bar=labor_bar)

        try
            sol = solve(model; init=init)
            init = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]

            if output == :real_gdp
                values[i] = real_gdp(sol)
            elseif output == :nominal_gdp
                values[i] = nominal_gdp(sol)
            elseif output == :sectoral_q
                sectoral_values[:, i] = sol.quantities
            elseif output == :sectoral_p
                sectoral_values[:, i] = sol.prices
            end
        catch e
            @warn "Solve failed at grid point $i (η=$η, ϵ=$ϵ, θ=$θ, σ=$σ): $e"
            values[i] = NaN
            if sectoral_values !== nothing
                sectoral_values[:, i] .= NaN
            end
        end
    end

    # ── Compute partial R² ──
    output_name = String(output)

    if output in (:real_gdp, :nominal_gdp)
        result = _compute_partial_r2(["η", "ϵ", "θ", "σ"], grid, values, output_name)
        return result
    else
        # Per-sector decomposition
        results = Vector{VarianceDecompositionResult}(undef, N)
        for s in 1:N
            results[s] = _compute_partial_r2(["η", "ϵ", "θ", "σ"], grid, sectoral_values[s, :], "$(output_name)_sector_$(s)")
        end
        return results
    end
end

"""
    _compute_partial_r2(factor_names, grid, values, output_name)

Internal helper: compute partial R² for each factor using ANOVA-style decomposition.
"""
function _compute_partial_r2(factor_names::Vector{String}, grid::DataFrame, values::Vector{Float64}, output_name::String)
    # Remove NaN values
    valid = .!isnan.(values)
    n_valid = sum(valid)

    if n_valid == 0
        return VarianceDecompositionResult(factor_names, Dict(f => NaN for f in factor_names), grid, values, output_name)
    end

    v = values[valid]
    grand_mean = mean(v)
    ss_total = sum((v .- grand_mean) .^ 2)

    if ss_total ≈ 0
        # No variation — all factors contribute nothing
        pr2 = Dict(f => 0.0 for f in factor_names)
        return VarianceDecompositionResult(factor_names, pr2, grid, values, output_name)
    end

    pr2 = Dict{String, Float64}()
    for f in factor_names
        col = grid[!, f][valid]
        levels = unique(col)
        ss_between = 0.0
        for lvl in levels
            mask = col .== lvl
            group_mean = mean(v[mask])
            n_group = sum(mask)
            ss_between += n_group * (group_mean - grand_mean)^2
        end
        pr2[f] = ss_between / ss_total
    end

    return VarianceDecompositionResult(factor_names, pr2, grid, values, output_name)
end

"""
    summary_table(vd::VarianceDecompositionResult)

Print a formatted summary of the variance decomposition.
"""
function summary_table(vd::VarianceDecompositionResult)
    println("\n═ Variance Decomposition: $(vd.output_name) ═")
    println("─" ^ 45)
    @printf("%-10s %12s %8s\n", "Factor", "Partial R²", "% Share")
    println("─" ^ 45)

    total_explained = sum(values(vd.partial_r2))
    for f in vd.factors
        val = vd.partial_r2[f]
        pct = total_explained > 0 ? 100 * val / total_explained : 0.0
        @printf("%-10s %12.4f %7.1f%%\n", f, val, pct)
    end
    println("─" ^ 45)
    @printf("%-10s %12.4f\n", "Total", total_explained)
    @printf("%-10s %12.4f\n", "Residual", 1.0 - total_explained)
    println("═" ^ 45)
end

"""
    eta_sweep_full(data, shocks; θ=0.5, ϵ=0.5, σ=0.9, labor_bar=nothing)

Run a comprehensive η sweep from perfectly inelastic to perfectly elastic labor.
Returns an `EtaSweepResult` with solutions at each η value.

This is the pilot sweep for the go/no-go decision.
"""
function eta_sweep_full(data::Data, shocks::Shocks; θ=0.5, ϵ=0.5, σ=0.9, labor_bar::Union{Float64, Nothing}=nothing)
    η_values = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0]
    sols = eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=labor_bar)
    return EtaSweepResult(η_values, sols)
end

"""
    pilot_eta_sweep(data, shocks; kwargs...)

Run the pilot η sweep and produce a diagnostic summary.

This implements the go/no-go decision criterion from the workplan:
- Check whether sectoral results show interesting variation across η
- Check whether η dominates (ε, θ, σ) in the variance decomposition

Prints diagnostics and returns both the sweep results and a recommendation.
"""
function pilot_eta_sweep(data::Data, shocks::Shocks; θ::Float64=0.5, ϵ::Float64=0.5, σ::Float64=0.9)
    println("\n" * "="^70)
    println("  PILOT η SWEEP — Go/No-Go Decision")
    println("="^70)

    # ── Part 1: η sweep ──
    println("\n[1/2] Running η sweep...")
    esr = eta_sweep_full(data, shocks; θ=θ, ϵ=ϵ, σ=σ)

    gdp_values = real_gdp_sweep(esr)
    println("\n  η values and real GDP:")
    for (i, η) in enumerate(esr.η_values)
        @printf("    η = %-8.2f  →  real GDP = %.6f  (Δ = %+.4f%%)\n",
            η, gdp_values[i], 100*(gdp_values[i] - gdp_values[1]))
    end

    # Check variation in sectoral quantities
    Q = sectoral_quantities(esr)
    sectoral_variation = vec(std(Q, dims=2))
    max_var = maximum(sectoral_variation)
    max_var_sector = argmax(sectoral_variation)
    sector_name = data.io.Sektoren[max_var_sector]

    println("\n  Sectoral variation across η:")
    @printf("    Max CV: sector %d (%s) = %.4f\n", max_var_sector, sector_name, max_var)
    @printf("    Mean CV: %.4f\n", mean(sectoral_variation))
    @printf("    Median CV: %.4f\n", median(sectoral_variation))

    # ── Part 2: Quick variance decomposition (reduced grid for speed) ──
    println("\n[2/2] Running variance decomposition (reduced grid)...")
    vd = variance_decomposition(data, shocks;
        η_values = [0.0, 0.5, 1.0, 5.0, 50.0],
        ϵ_values = [0.1, 0.5, 0.99],
        θ_values = [0.5],  # reduced for speed
        σ_values = [0.5],  # reduced for speed
        output = :real_gdp,
        verbose = false,
    )
    summary_table(vd)

    # ── Go/No-Go assessment ──
    η_share = vd.partial_r2["η"]
    other_share = sum([vd.partial_r2["ϵ"], vd.partial_r2["θ"], vd.partial_r2["σ"]])

    println("\n" * "="^70)
    println("  ASSESSMENT")
    println("="^70)

    interesting_sectors = count(v -> v > 0.01 * max_var, sectoral_variation)
    println("  Sectors with >1% of max variation: $interesting_sectors / $(length(sectoral_variation))")

    if η_share > 0
        @printf("  η accounts for %.1f%% of explained variance\n", 100 * η_share / (η_share + other_share))
    end

    go = (interesting_sectors >= 5) && (max_var > 0.01) && (η_share > 0.05)
    println("\n  Recommendation: $(go ? "GO ✅" : "NO-GO ⚠️")")
    if !go
        println("    Reason: ")
        if interesting_sectors < 5
            println("    - Too few sectors show interesting variation ($interesting_sectors < 5)")
        end
        if max_var <= 0.01
            println("    - Maximum sectoral variation too small ($(max_var))")
        end
        if η_share <= 0.05
            println("    - η contribution to variance too small ($(η_share))")
        end
    end

    println("="^70)

    return esr, vd, go
end
