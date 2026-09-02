# ═══════════════════════════════════════════════════════════════════════════════
# Variance Decomposition — Sobol Indices on a Full Factorial Grid
# ═══════════════════════════════════════════════════════════════════════════════
#
# Implements:
#   1. η sweep — run the mobile-labor CES model across a vector of η values
#   2. Sobol variance decomposition — factorial design over (η, ε, θ, σ) with
#      first-order and total-order Sobol indices
#
# Method:
#   Sobol (1993) decomposition on a balanced full factorial grid.
#   - First-order index S_f = SS_f / SS_total
#     (variance attributable to factor f alone)
#   - Total-order index ST_f = 1 − SS_{-f} / SS_total
#     (variance that requires factor f in any form — main effect or interaction)
#   - Interaction share = ST_f − S_f
#     (variance only explained by f interacting with other factors)
#   - Unexplained = 1 − sum(S_f)  (should be zero in a full factorial;
#     any positive value indicates missing grid points)
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

using ProgressMeter
using Printf
using Statistics
using CSV
using DataFrames

"""
    eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=nothing)

Run the mobile-labor CES model across a vector of η values.
Returns a vector of `Solution` objects.

# Example
```julia
η_grid = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
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
    SobolResult

Results of a Sobol variance decomposition on a full factorial grid.

- `factors`: names of the factors (e.g., ["η", "ϵ", "θ", "σ"])
- `S_f`: first-order Sobol index for each factor (main effect, S_f = SS_f / SS_total)
- `ST_f`: total-order Sobol index for each factor (ST_f = 1 − SS_{-f} / SS_total)
- `grid`: the full factorial grid used
- `values`: the output values at each grid point (NaN where solver failed)
- `n_failed`: number of grid points with NaN (solver failure)
- `output_name`: name of the output variable decomposed
- `ss_total`: total sum of squares
- `ss_explained`: sum of main-effect SS (sum(S_f) × SS_total)
- `frac_unexplained`: 1 − sum(S_f) (= interaction share in balanced data)
"""
struct SobolResult
    factors::Vector{String}
    S_f::Dict{String, Float64}
    ST_f::Dict{String, Float64}
    grid::DataFrame
    values::Vector{Float64}
    n_failed::Int
    output_name::String
    ss_total::Float64
    ss_explained::Float64
    frac_unexplained::Float64
end

"Deprecated compatibility alias. Use `SobolResult`."
const VarianceDecompositionResult = SobolResult

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
- `SobolResult` with first-order and total-order indices for each factor

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
    y_vals = Vector{Float64}(undef, n)
    sectoral_y_vals = output in (:sectoral_q, :sectoral_p) ? zeros(N, n) : nothing

    init = [ones(N); data.λ; 1.0]

    @showprogress "Variance decomposition: " for i in 1:n
        η = grid.η[i]
        ϵ = grid.ϵ[i]
        θ = grid.θ[i]
        σ = grid.σ[i]

        # Preserve the historical finite approximation for an infinite
        # extrapolation request. This is not an elastic-labor-supply limit.
        if isinf(η)
            # Approximate the geometric reallocation exponent with a large value.
            η_eff = 1e6
        else
            η_eff = η
        end

        model = mobile_labor_model(data, shocks, θ, ϵ, σ, η_eff; labor_bar=labor_bar)

        try
            sol = solve(model; init=init)
            init = [sol.prices_raw; sol.quantities; sol.wages_raw[1]]

            if output == :real_gdp
                y_vals[i] = real_gdp(sol)
            elseif output == :nominal_gdp
                y_vals[i] = nominal_gdp(sol)
            elseif output == :sectoral_q
                sectoral_y_vals[:, i] = sol.quantities
            elseif output == :sectoral_p
                sectoral_y_vals[:, i] = sol.prices
            end
        catch e
            @warn "Solve failed at grid point $i (η=$η, ϵ=$ϵ, θ=$θ, σ=$σ): $e"
            y_vals[i] = NaN
            if sectoral_y_vals !== nothing
                sectoral_y_vals[:, i] .= NaN
            end
        end
    end

    # ── Compute Sobol indices ──
    output_name = String(output)

    if output in (:real_gdp, :nominal_gdp)
        result = _compute_sobol_indices(["η", "ϵ", "θ", "σ"], grid, y_vals, output_name)
        return result
    else
        # Per-sector decomposition (kept for backward compatibility)
        results = Vector{SobolResult}(undef, N)
        for s in 1:N
            results[s] = _compute_sobol_indices(["η", "ϵ", "θ", "σ"], grid, sectoral_y_vals[s, :], "$(output_name)_sector_$(s)")
        end
        return results
    end
end


"""
    _compute_sobol_indices(factor_names, grid, values, output_name)

Compute first-order (S_f) and total-order (ST_f) Sobol indices on a
balanced full factorial grid using ANOVA sum-of-squares decomposition.

First-order:   S_f = SS_f / SS_total
Total-order:   ST_f = 1 − SS_{-f} / SS_total

where SS_{-f} is the sum of squares explained by all factors EXCEPT f.

NaN values (from solver failures) are dropped; `n_failed` records how many
points were lost. The total-order index is biased if grid points are missing,
so `n_failed` should be reported.
"""
function _compute_sobol_indices(factor_names::Vector{String}, grid::DataFrame, output_vals::Vector{Float64}, output_name::String)
    n_total = length(output_vals)
    valid = .!isnan.(output_vals)
    n_valid = sum(valid)
    n_failed = n_total - n_valid

    if n_valid == 0
        S_f = Dict(f => NaN for f in factor_names)
        ST_f = Dict(f => NaN for f in factor_names)
        return SobolResult(factor_names, S_f, ST_f, grid, output_vals, n_failed, output_name, 0.0, 0.0, NaN)
            end

            v = output_vals[valid]
            grand_mean = mean(v)
            ss_total = sum((v .- grand_mean) .^ 2)

            if ss_total ≈ 0
                S_f = Dict(f => 0.0 for f in factor_names)
                ST_f = Dict(f => 0.0 for f in factor_names)
                return SobolResult(factor_names, S_f, ST_f, grid, output_vals, n_failed, output_name, ss_total, 0.0, NaN)
            end

            # ── First-order indices S_f ──
            S_f = Dict{String, Float64}()
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
                S_f[f] = ss_between / ss_total
            end

            # ── Total-order indices ST_f = 1 − SS_{-f} / SS_total ──
            # SS_{-f} = variance explained by all factors EXCEPT f.
            # On a balanced full factorial, for each combination of levels of all
            # other factors, we take the mean of y across f's levels, and sum the
            # squared deviations weighted by the number of f levels.
            ST_f = Dict{String, Float64}()
            for f in factor_names
                other = setdiff(factor_names, [f])
                # Group by all other factors
                other_cols = [grid[!, f_name][valid] for f_name in other]
                # Unique level combinations of the other factors
                combos = unique([Tuple([c[i] for c in other_cols]) for i in 1:n_valid])
                n_levels_f = length(unique(grid[!, f][valid]))
                ss_except_f = 0.0
                for combo in combos
                    mask = trues(n_valid)
                    for (j, f_name) in enumerate(other)
                        mask = mask .& (other_cols[j] .== combo[j])
                    end
                    n_in_mask = sum(mask)
                    if n_in_mask > 0
                        group_mean = mean(v[mask])
                        # Weight by number of f levels (what we'd expect in a balanced design)
                        ss_except_f += n_levels_f * (group_mean - grand_mean)^2
                    end
                end
                ST_f[f] = 1.0 - ss_except_f / ss_total
            end

            # Clip to [0, 1] for numerical safety
            for f in factor_names
                S_f[f] = max(0.0, min(1.0, S_f[f]))
                ST_f[f] = max(0.0, min(1.0, ST_f[f]))
            end

            ss_explained = sum(Base.values(S_f)) * ss_total
            frac_unexplained = 1.0 - sum(Base.values(S_f))

            return SobolResult(factor_names, S_f, ST_f, grid, output_vals, n_failed, output_name, ss_total, ss_explained, frac_unexplained)
end


"""
    summary_table(vd::SobolResult; save_csv=nothing)

Print a formatted summary of the Sobol variance decomposition.
Reports first-order (S_f) and total-order (ST_f) indices as absolute
shares (no renormalization). Also reports the interaction strength
ST_f − S_f and the unexplained fraction.

If `save_csv` is a file path, writes the results to CSV.
"""
function summary_table(vd::SobolResult; save_csv=nothing)
    println("\n═ Sobol Variance Decomposition: $(vd.output_name) ═")
    println("Grid: $(nrow(vd.grid)) points, $(vd.n_failed) failed (NaN) = $(round(100*vd.n_failed/nrow(vd.grid), digits=1))%")
    if vd.n_failed > 0
        println("⚠  Total-order indices are biased by missing grid points — interpret with caution")
    end
    println("─" ^ 62)
    @printf("%-10s %12s %12s %12s\n", "Factor", "S_f", "ST_f", "ST_f−S_f")
    @printf("         %12s %12s %12s\n", "(first-order)", "(total-order)", "(interaction)")
    println("─" ^ 62)

    for f in vd.factors
        sf = vd.S_f[f]
        st = vd.ST_f[f]
        inter = st - sf
        @printf("%-10s %12.4f %12.4f %12.4f\n", f, sf, st, inter)
    end
    println("─" ^ 62)
    @printf("%-10s %12.4f\n", "Sum S_f", sum(values(vd.S_f)))
    @printf("%-10s %12.4f\n", "1−Sum S_f", vd.frac_unexplained)
    @printf("%-10s %12.4f\n", "SS_total", vd.ss_total)
    println("═" ^ 62)

    # Interpretation guide
    if vd.n_failed == 0
        println("Interpretation: 1−Sum S_f = interaction share (balanced factorial).")
    else
        println("Interpretation: 1−Sum S_f partly includes missing-data bias.")
    end
    dom = argmax([vd.S_f[f] for f in vd.factors])
    @printf("Dominant factor: %s (S_f=%.4f)\n", vd.factors[dom], vd.S_f[vd.factors[dom]])

    # ── Save CSV ──
    if save_csv !== nothing
        df = DataFrame(
            Factor = vd.factors,
            S_f = [vd.S_f[f] for f in vd.factors],
            ST_f = [vd.ST_f[f] for f in vd.factors],
            Interaction = [vd.ST_f[f] - vd.S_f[f] for f in vd.factors],
        )
        push!(df, ("Sum", sum(values(vd.S_f)), sum(values(vd.ST_f)), sum(values(vd.ST_f)) - sum(values(vd.S_f))))
        push!(df, ("1−Sum_S_f", vd.frac_unexplained, NaN, NaN))
        push!(df, ("SS_total", vd.ss_total, NaN, NaN))
        push!(df, ("n_failed", Float64(vd.n_failed), NaN, NaN))
        push!(df, ("n_grid", Float64(nrow(vd.grid)), NaN, NaN))
        CSV.write(save_csv, df)
        @printf("Results saved to %s\n", save_csv)
    end
    println()
end

"""
    eta_sweep_full(data, shocks; θ=0.5, ϵ=0.5, σ=0.9, labor_bar=nothing)

Run a comprehensive sweep over the intersectoral reallocation parameter η.
Returns an `EtaSweepResult` with solutions at each η value.

This is a descriptive sweep; it does not make a go/no-go claim.
"""
function eta_sweep_full(data::Data, shocks::Shocks; θ=0.5, ϵ=0.5, σ=0.9, labor_bar::Union{Float64, Nothing}=nothing)
    η_values = [0.0, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0]
    sols = eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=labor_bar)
    return EtaSweepResult(η_values, sols)
end

"""Run the sweep and return descriptive diagnostics, without a binary decision."""
function eta_sweep_diagnostics(data::Data, shocks::Shocks; θ::Float64=0.5, ϵ::Float64=0.5, σ::Float64=0.9, labor_bar::Union{Float64, Nothing}=nothing)
    esr = eta_sweep_full(data, shocks; θ=θ, ϵ=ϵ, σ=σ, labor_bar=labor_bar)
    variation = vec(std(sectoral_quantities(esr), dims=2))
    vd = variance_decomposition(data, shocks; η_values=[0.0, 0.5, 1.0, 5.0, 50.0],
        ϵ_values=[0.1, 0.5, 0.99], θ_values=[0.5], σ_values=[0.5], output=:real_gdp,
        labor_bar=labor_bar, verbose=false)
    (sweep=esr, decomposition=vd, sectoral_variation=variation,
     max_variation=maximum(variation), eta_share=vd.S_f["η"],
     other_share=sum(vd.S_f[f] for f in ("ϵ", "θ", "σ")))
end

"""
    pilot_eta_sweep(data, shocks; kwargs...)

Deprecated compatibility wrapper for the former pilot/go-no-go report.
"""
function pilot_eta_sweep(data::Data, shocks::Shocks; θ::Float64=0.5, ϵ::Float64=0.5, σ::Float64=0.9)
    Base.depwarn("pilot_eta_sweep is deprecated; use eta_sweep_diagnostics", :pilot_eta_sweep)
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
    η_share = vd.S_f["η"]
    other_share = sum([vd.S_f["ϵ"], vd.S_f["θ"], vd.S_f["σ"]])

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
