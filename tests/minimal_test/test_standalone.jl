# ═══════════════════════════════════════════════════════════════════════════════
# Self-contained test — uses ONLY Julia stdlib, no Pkg.instantiate needed
# ═══════════════════════════════════════════════════════════════════════════════
# Implements a simple Newton solver to avoid the NonlinearSolve dependency.
# Tests the economic logic of the mobile labor model.

using LinearAlgebra
using Printf
using Statistics
using DelimitedFiles

cd("/workspace/BFrep/(3)BeyondHulten")

# ═══════════════════════════════════════════════════════════════════════════════
# Data loading (minimal reimplementation of read_data/generate_data)
# ═══════════════════════════════════════════════════════════════════════════════

println("Loading data...")

# Read CSV the hard way (no CSV.jl)
# The file has a 2-line header: line 1 = "Jahr (starts a quoted field)
# line 2 = Gütergruppen, ...";col2;col3;... (completes the quoted field)
function read_io_table(filename)
    raw = read(filename, String)
    # Normalize line endings
    raw = replace(raw, "\r\n" => "\n")
    lines = split(raw, "\n")

    # Header is lines 1+2 combined
    header_line = lines[1] * "\n" * lines[2]
    header = split(header_line, ';')
    # First header cell is the quoted "Jahr\n..." — strip quotes
    header[1] = replace(header[1], "\"" => "")
    n_cols = length(header)

    # Data rows start from line 3
    data_lines = filter(!isempty, lines[3:end])

    sectors = String[]
    data_matrix = zeros(length(data_lines), n_cols - 1)

    for (i, row) in enumerate(data_lines)
        parts = split(row, ';')
        push!(sectors, strip(parts[1]))
        for j in 2:min(n_cols, length(parts))
            val_str = replace(parts[j], "," => ".")
            val_str = strip(val_str)
            if val_str in ["-", "x", ""]
                data_matrix[i, j-1] = 0.0
            else
                data_matrix[i, j-1] = parse(Float64, val_str)
            end
        end
    end

    return sectors, header[2:end], data_matrix
end

sectors, colnames, io_matrix = read_io_table("data/I-O_DE2019_formatiert.csv")
println("  IO table: $(length(sectors)) rows × $(length(colnames)) cols")
println("  First 5 sectors: $(sectors[1:5])")

# Strip whitespace from colnames for matching
colnames = strip.(colnames)
sectors = strip.(sectors)

# Find key rows (sector names in first column)
N = 71
bws_row = findfirst(==("Bruttowertschöpfung"), sectors)
prodwert_row = findfirst(==("Produktionswert"), sectors)
arbeit_row = findfirst(==("Arbeitnehmerentgelt im Inland"), sectors)

# Find key columns
gesamt_verwendung_col = findfirst(==("Gesamte Verwendung von Gütern"), colnames)
konsum_col = findfirst(==("Konsumausgaben der privaten Haushalte im Inland"), colnames)
exporte_col = findfirst(==("Exporte"), colnames)
letzte_verwendung_col = findfirst(==("Letzte Verwendung von Gütern zusammen"), colnames)

println("  Bruttowertschöpfung row: $bws_row")
println("  Produktionswert row: $prodwert_row")
println("  Arbeitnehmerentgelt row: $arbeit_row")
println("  Gesamte Verwendung col: $gesamt_verwendung_col")
println("  Konsum col: $konsum_col")
println("  Exporte col: $exporte_col")

# Extract model data (mirrors generate_data)
Ω = io_matrix[1:N, 1:N]
Ω = Ω ./ sum(Ω, dims=2)
Ω[isnan.(Ω)] .= 0.0

grossy = io_matrix[1:N, gesamt_verwendung_col]
value_added = io_matrix[bws_row, 1:N]
labor_comp = io_matrix[arbeit_row, 1:N]

factor_share = value_added ./ grossy
# Clamp to valid range
factor_share = clamp.(factor_share, 0.01, 0.99)

# Consumption: sum across consumption categories
consumption = zeros(N)
for col in konsum_col:exporte_col
    consumption .+= io_matrix[1:N, col]
end

consumption_share_gross_output = consumption ./ grossy

# Compute consumption_share via Leontief inverse
consumption_share_raw = (I - diagm(1 .- factor_share) * Ω)' * grossy
consumption_share_raw[consumption_share_raw .< 0] .= 0.0
consumption_share = consumption_share_raw ./ sum(consumption_share_raw)

# λ = Leontief inverse × consumption_share
λ = inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share
labor_share = λ .* factor_share

println("\n  factor_share[1:5]: $(round.(factor_share[1:5], digits=3))")
println("  sum(consumption_share): $(sum(consumption_share))")
println("  sum(λ): $(sum(λ))")
println("  sum(labor_share): $(sum(labor_share))")
println("  grossy[1:5]: $(round.(grossy[1:5], digits=1))")

# ═══════════════════════════════════════════════════════════════════════════════
# Standard shock
# ═══════════════════════════════════════════════════════════════════════════════
target_sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"
sector_idx = findfirst(==(target_sector), sectors[1:N])
demand_shock = ones(N)
supply_shock = ones(N)
demand_shock[sector_idx] = 1.8097957577943152
demand_shock_raw = zeros(N)
demand_shock_raw[sector_idx] = 100_000.0  # 100k thousand € = 100M €

println("\n  Shock sector: $target_sector (idx=$sector_idx)")
println("  demand_shock[sector_idx]: $(demand_shock[sector_idx])")

# ═══════════════════════════════════════════════════════════════════════════════
# Mobile labor CES model (inline implementation)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    mobile_labor_problem(X, params)

Evaluate the residual of the mobile-labor CES equilibrium system.

Unknowns: X = [p(1:N); y(1:N); w]  — 2N+1 elements
Equations (2N+1):
  1. Price equations        (N-1):  p_i = cost(p, w, P_i)   for i=2..N
  2. Market clearing        (N):     y_i = D_intermediary + D_final
  3. Labor market clearing  (1):     Σ L_i(p,y,w) = L̄ · w^η
  4. Numeraire              (1):     Σ β_i · p_i^(1-σ) = 1   (CPI = 1)

Note: We drop the first price equation (i=1) and replace it with the
numeraire constraint. This breaks the price-level indeterminacy that
would otherwise make the wage unresponsive to η.
"""
function mobile_labor_problem(X, params)
    (; Ω, consumption_share, factor_share, supply_shock, demand_shock, labor_bar, θ, ϵ, σ, η) = params
    N = length(factor_share)

    p = max.(X[1:N], 1e-10)
    y = max.(X[N+1:2N], 1e-10)
    w = max(X[2N+1], 1e-10)

    out = zeros(eltype(X), 2N + 1)

    # Intermediate goods price index
    intermediate_price = (Ω * p .^ (1 - θ)) .^ (1 / (1 - θ))

    # CPI (for numeraire and final demand)
    cpi = sum(consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))

    # Sectoral labor demand: L_i = [p_i · A_i^((ε-1)/ε) · α_i^(1/ε) · y_i^(1/ε) / w]^ε
    L_i = (p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) ./ w) .^ ϵ

    # Final demand
    total_income = w * sum(L_i)
    final_demand = (total_income * p .^ (-σ) .* demand_shock .* consumption_share) ./ cpi .^ (-σ)

    # Intermediary demand
    intermediary_demand = p .^ (-θ) .* (Ω' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))

    # Eq 1: Price equations for sectors 2..N (drop sector 1 → numeraire)
    cost = (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) .+ (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
    out[1:N-1] = p[2:N] .- cost[2:N]

    # Eq 2: Market clearing (N equations)
    out[N:2N-1] = y .- intermediary_demand .- final_demand

    # Eq 3: Labor market clearing
    out[2N] = sum(L_i) - labor_bar * w^η

    # Eq 4: Numeraire constraint — CPI = 1
    out[2N+1] = cpi - 1.0

    return out
end

"""
    numerical_jacobian(f, x, params; ε=1e-6)

Compute the Jacobian of f(x, params) via finite differences.
"""
function numerical_jacobian(f, x, params; eps=1e-6)
    n = length(x)
    f0 = f(x, params)
    m = length(f0)
    J = zeros(m, n)
    for j in 1:n
        x_plus = copy(x)
        x_plus[j] += eps
        f_plus = f(x_plus, params)
        J[:, j] = (f_plus - f0) / eps
    end
    return J
end

"""
    solve_newton(f, x0, params; max_iter=100, tol=1e-7, verbose=false)

Simple damped Newton solver with fallback to fixed-point iteration.
"""
function solve_newton(f, x0, params; max_iter=200, tol=1e-7, verbose=false)
    x = copy(x0)
    n = length(x)

    for iter in 1:max_iter
        F = f(x, params)
        norm_F = norm(F)

        if verbose && (iter % 10 == 0 || iter == 1)
            @printf("    iter %d: ||F|| = %.2e\n", iter, norm_F)
        end

        if norm_F < tol
            if verbose
                @printf("    Converged at iter %d: ||F|| = %.2e\n", iter, norm_F)
            end
            return x, true
        end

        # Jacobian
        J = numerical_jacobian(f, x, params)

        # Newton step with damping
        try
            dx = -(J \ F)
            # Damping: don't take too large steps
            alpha = min(1.0, 1.0 / max(1.0, norm(dx) / norm(x)))
            x_new = x + alpha * dx

            # Check if step improved
            F_new = f(x_new, params)
            if norm(F_new) < norm_F
                x = x_new
            else
                # Smaller step
                x_new = x + 0.1 * alpha * dx
                x = x_new
            end
        catch e
            # Singular Jacobian — perturb
            if verbose
                println("    Jacobian singular, perturbing...")
            end
            x .+= 1e-4 * randn(n)
        end
    end

    return x, false
end

# ═══════════════════════════════════════════════════════════════════════════════
# Run tests
# ═══════════════════════════════════════════════════════════════════════════════

labor_bar = sum(labor_share)

println("\n" * "="^70)
println("  Testing Mobile Labor CES Model")
println("="^70)

# ── Test across η values ──
η_values = [0.0, 0.01, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
results = []

# Common parameters
θ_val = 0.5
ϵ_val = 0.5
σ_val = 0.9

params_base = (
    Ω = Ω,
    consumption_share = consumption_share,
    factor_share = factor_share,
    supply_shock = supply_shock,
    demand_shock = demand_shock,
    labor_bar = labor_bar,
    θ = θ_val,
    ϵ = ϵ_val,
    σ = σ_val,
)

println("\nRunning η sweep:")
println("  η        Real GDP    Wage      ||F||      Converged")
println("  " * "-"^60)

# Initial guess
x0 = [ones(N); λ; 1.0]

for η in η_values
    params = merge(params_base, (η = η,))
    x_init = copy(x0)

    x_sol, converged = solve_newton(mobile_labor_problem, x_init, params, max_iter=300, tol=1e-6)

    p_sol = x_sol[1:N]
    q_sol = x_sol[N+1:2N]
    w_sol = x_sol[2N+1]

    # Compute real GDP (laspeyres index)
    consumption_adj = demand_shock .* consumption_share
    numeraire = (consumption_share' * p_sol .^ (1 - σ_val))^(1 / (1 - σ_val))
    L_i_sol = (p_sol .* (supply_shock .^ ((ϵ_val - 1) / ϵ_val)) .* (factor_share .^ (1 / ϵ_val)) .* (q_sol .^ (1 / ϵ_val)) ./ w_sol) .^ ϵ_val
    total_income = w_sol * sum(L_i_sol)
    consumption_sol = total_income .* consumption_adj .* (p_sol / numeraire) .^ (-σ_val)
    laspeyres = sum(consumption_sol) / sum(consumption_share)
    nominal_gdp = total_income / numeraire

    push!(results, (η, laspeyres, nominal_gdp, w_sol, converged, p_sol, q_sol))

    @printf("  %-8.2f %.6f   %.6f  %s  %s\n",
        η, laspeyres, w_sol,
        @sprintf("%.2e", norm(mobile_labor_problem(x_sol, params))),
        converged ? "✅" : "⚠️")

    # Warm start next
    global x0 = [p_sol; q_sol; w_sol]
end

# ── Sectoral variation analysis ──
println("\n" * "="^70)
println("  Sectoral Variation Analysis")
println("="^70)

# Extract quantity matrix (sectors × η)
Q_mat = hcat([r[7] for r in results]...)
cv_sectors = vec(std(Q_mat, dims=2))

@printf("\n  Max CV:  sector %d (%s) = %.4f\n", argmax(cv_sectors), sectors[argmax(cv_sectors)], maximum(cv_sectors))
@printf("  Mean CV: %.4f\n", mean(cv_sectors))
@printf("  Median CV: %.4f\n", median(cv_sectors))
@printf("  # sectors with CV > 0.001: %d / %d\n", count(cv_sectors .> 0.001), N)
@printf("  # sectors with CV > 0.01:  %d / %d\n", count(cv_sectors .> 0.01), N)

# Top 5 most varying sectors
sorted_idx = sortperm(cv_sectors, rev=true)
println("\n  Top 5 most varying sectors:")
for i in sorted_idx[1:5]
    @printf("    %-50s CV = %.4f\n", sectors[i], cv_sectors[i])
end

# ── GDP range analysis ──
println("\n" * "="^70)
println("  GDP Range Analysis")
println("="^70)

gdp_values = [r[2] for r in results]
gdp_min, gdp_max = extrema(gdp_values)
@printf("  GDP range: [%.6f, %.6f], spread = %.4f%%\n", gdp_min, gdp_max, 100*(gdp_max - gdp_min))
@printf("  GDP at η=0:   %.6f\n", gdp_values[1])
@printf("  GDP at η=1:   %.6f\n", findfirst(==(1.0), η_values) |> i -> gdp_values[i])
@printf("  GDP at η=∞:   %.6f (η=100)\n", gdp_values[end])

println("\n" * "="^70)
println("  Done!")
println("="^70)