# src/core/diagnostics.jl — utilities, impulses, variance decomposition
# (concatenated verbatim: src/util.jl + src/impulses.jl + src/variance_decomposition.jl, in that order)

using ProgressMeter, ThreadsX


"""Construct the standard demand shock for one modeled sector."""
function standard_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	demand_shock = ones(n)
	supply_shock = ones(n)
	demand_shock[index] = 1.8097957577943152
	shocks = Shocks(supply_shock, demand_shock, zeros(n))
	return shocks
end


"""
	autonomous_shock(data; sector, autonomous_mult, investment_mult)

Construct additive autonomous and investment-demand multipliers for `sector`.
Inside `MobileLaborCES`, each multiplier is scaled by that sector's consumption
share and baseline aggregate labor income. Other sectors receive zero additive
demand.
"""
function autonomous_shock(data;
	sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten",
	autonomous_mult = 1.8097957577943152,
	investment_mult = 0.0)
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	aut = zeros(n)
	aut[index] = autonomous_mult
	inv = zeros(n)
	inv[index] = investment_mult
	Shocks(ones(n), ones(n), zeros(n), aut, inv)
end


"""Construct the standard technology shock for one modeled sector."""
function standard_tech_shock(data, sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten")
	n = length(data.grossy)
	index = findfirst(==(sector), data.io.Sektoren)
	(index === nothing || index > n) && throw(ArgumentError("sector $sector is not a modeled sector"))
	demand_shock = ones(n)
	supply_shock = ones(n)
	supply_shock[index] = 1.2
	Shocks(supply_shock, demand_shock, zeros(n))
end

"""Construct demand and supply shocks from an impulse-response table."""
function impulse_shock(data, impulses)
	# `impulses` contains year, one column per modeled sector, and wages. Retain every goods
	# sector while excluding the non-sector columns at either end.
	impules_2019_prices = impulses[:, 2:end-1] ./ inflator
	size(impules_2019_prices, 2) == length(data.grossy) ||
		throw(DimensionMismatch("the impulse table must contain one column per modeled sector"))
	n = length(data.grossy)
	last_use = Vector(data.io[1:n, "Letzte Verwendung von Gütern zusammen"])
	effect = 1 .+ impules_2019_prices ./ last_use'
	demand_shock = [mean(col) for col in eachcol(effect[1:min(2, size(effect, 1)), :])]
	supply_shock = ones(n)
	Shocks(supply_shock, demand_shock, [mean(col) for col in eachcol(impules_2019_prices)])
end
struct ElasticityGradientSolution
	ϵ::Vector{Solution}
	θ::Vector{Solution}
	σ::Vector{Solution}
	labor_realloc::Bool
	nominal::Bool
end


function gradient(data, shocks, labor_slack, labor_reallocation, elasticity, sol, el, nominal = false)::Vector{Solution}
	len = 1000
	sols = Vector{Solution}(undef, len)
	arr = copy(el)
	u0 = [sol.prices; sol.quantities]
	@inbounds for (idx, i) in enumerate(range(0.99, 0.015, len))
		arr[elasticity] = i
		elasticities = CESElasticities(arr...)
		ces = CES(elasticities, labor_slack, labor_reallocation)
		model = Model(data, shocks, ces)
		sol_prev = solve(model, init = u0)
		u0 = [sol_prev.prices_raw; sol_prev.quantities]
		sols[idx] = sol_prev
	end
	return sols
end

function elasticity_gradient(data,
	shocks,
	labor_slack = full_labor_slack,
	labor_reallocation = false,
	starting_elasticities = [0.99, 0.99, 0.99],
	nominal = false,
)


	elasticities = CESElasticities(starting_elasticities...)
	ces = CES(elasticities, labor_slack, labor_reallocation)
	model = Model(data, shocks, ces)
	sol_original = solve(model)
	sols_ϵ, sols_θ, sols_σ = fetch.(Threads.@spawn(gradient($data, $shocks, $labor_slack, $labor_reallocation, i, $sol_original, $starting_elasticities, $nominal)) for i in 1:3)

	return ElasticityGradientSolution(sols_ϵ, sols_θ, sols_σ, labor_reallocation, nominal)
end


"""
Returns the consumer price index (β ̇ p^(1 - σ))^(1/(1-σ) 

"""
function cpi(sol::Solution)
	σ = sol.model.options.elasticities.σ
	(sol.model.data.consumption_share' * sol.prices_raw .^ (1- σ))^(1/(1 - σ))
end

# --- National-accounts GDP measurement (ADR-0018) ---

"""
	gdp_components(model::Model{MobileLaborCES}, sol::Solution) -> NamedTuple

Seven signed aggregate GDP components at a solved open-economy equilibrium,
with `p = sol.prices_raw`, `y = sol.quantities`, `w = sol.wages_raw[1]`,
`F = sol.external_transfer` and
the demand blocks of the shared `_mobile_market_demand` hook
(`c_dom`, `additive`, `L_i`, `E`):

1. `C_gross`: gross household consumption `c_dom ./ (1 .- m)`;
2. `G+programme`: government demand plus the financed programme bundle;
3. `I`: exogenous investment; 4. `X`: exports (no margin);
5. `−M_final`: final-import margin content of C+G+I (including the programme);
6. `−M_int`: intermediate imports (row 74, ADR-0012);
7. `−T_int`: product taxes on intermediate use (row 75, ADR-0013).

The intermediate-bill leaks (6-7) are valued with the ADR-0016 CES bill
factor `k = p^ϵ · a^(ϵ−1) · P^(1−ϵ)` (same factor as
`external_balance_canary`, duplicated here so the canary is untouched).

Returns `(V, Q, wage_bill, wedge, external_transfer, programme_financing,
external_financing)` with `V`/`Q` the seven component values / quantities
(model units), `wage_bill = w·ΣL_i` and `wedge = sum(V) - wage_bill`. The
demand hook is evaluated with the solution's external transfer
(`external_transfer = sol.external_transfer`); `programme_financing =
Σ p_i g_i` under ExternalDebt (F3) and 0 otherwise, and `external_financing
= external_transfer + programme_financing` is the booked external position.
The exact accounting identity (ADR-0019) is `wedge = −canary.diff`: the
wedge IS the negated external-account identity gap, on and off equilibrium.
It is ≈ 0 at every η = 1 solution (mobile and fixed-wage, where the
cost-minimizing allocation lets zero-profit plus clearing close the external
account); at BF η = 0 it carries the fixed-allocation factor-market gap
(zero-profit prices the cost-minimizing labour demand, not the frozen
baseline allocation). Throws `ArgumentError` for other model types (the
closed cores have no open-economy blocks).
"""
function gdp_components(model::Model{MobileLaborCES}, sol::Solution)
	p = sol.prices_raw
	y = sol.quantities
	w = sol.wages_raw[1]
	blocks = _mobile_market_demand(model, p, y, w; external_transfer = sol.external_transfer)
	(; data, options, shocks) = model
	(; θ, ϵ) = options.elasticities
	m = data.import_margin
	gov = data.gov_demand
	add = blocks.additive
	exo = data.exo_demand
	x = data.exports_demand
	# 1. Gross household consumption (domestic block grossed up by the margin).
	c_gross = blocks.c_dom ./ (1 .- m)
	V1 = dot(p, c_gross)
	Q1 = sum(c_gross)
	# 2. Government consumption plus the financed programme bundle.
	V2 = dot(p, gov .+ add)
	Q2 = sum(gov .+ add)
	# 3. Investment. 4. Exports (domestic sales abroad: no margin).
	V3 = dot(p, exo)
	Q3 = sum(exo)
	V4 = dot(p, x)
	Q4 = sum(x)
	# 5. Final-import margin content of household and injected demand.
	M_cons = dot(p .* (m ./ max.(1 .- m, eps(Float64))), blocks.c_dom)
	M_inj = dot(p .* m, add .+ gov .+ exo)
	V5 = -(M_cons + M_inj)
	Q5 = sum(m .* c_gross) + sum(m .* (add .+ gov .+ exo))
	# 6-7. Intermediate-bill leaks, valued with the ADR-0016 CES bill factor
	# (duplicated from `external_balance_canary` so the canary is untouched).
	k = p .^ ϵ .* shocks.supply_shock .^ (ϵ - 1) .* _intermediate_price(data.Ω_raw, p, θ) .^ (1 - ϵ)
	V6 = -dot(k .* (data.M_int ./ data.λ), y)
	Q6 = sum((data.M_int ./ data.λ) .* y)
	V7 = -dot(k .* (data.T_int ./ data.λ), y)
	Q7 = sum((data.T_int ./ data.λ) .* y)
	V = Float64[V1, V2, V3, V4, V5, V6, V7]
	Q = Float64[Q1, Q2, Q3, Q4, Q5, Q6, Q7]
	wage_bill = w * sum(blocks.L_i)
	external_transfer = Float64(sol.external_transfer)
	programme_financing = model.financing isa ExternalDebt ? dot(p, blocks.additive) : 0.0
	external_financing = external_transfer + programme_financing
	return (; V = V, Q = Q, wage_bill = Float64(wage_bill), wedge = sum(V) - wage_bill,
		external_transfer = external_transfer, programme_financing = programme_financing,
		external_financing = external_financing)
end

function gdp_components(model::Model, ::Solution)
	throw(ArgumentError("gdp_components is only defined for the open-economy Model{MobileLaborCES}; " *
		"got $(typeof(model.options)) (the closed cores have no open-economy blocks)"))
end

"""
	gdp_deflator(sol::Solution, base::Solution) -> Float64

Törnqvist price index of the seven `gdp_components` aggregates (ADR-0018):
unit values `|V_j|/Q_j` with signed value shares `V_j/ΣV`, 1.0 at `base`.
A component whose value is zero at either end contributes no measured price
change; components that are zero at both ends are dropped; slots with
non-positive quantities are skipped (never divided by zero). Throws
`DomainError` if either side's total signed value is non-positive.
"""
function gdp_deflator(sol::Solution, base::Solution)::Float64
	c = gdp_components(sol.model, sol)
	c0 = gdp_components(base.model, base)
	tot = sum(c.V)
	tot0 = sum(c0.V)
	tot0 <= 0 && throw(DomainError(tot0, "base GDP value must be positive"))
	tot <= 0 && throw(DomainError(tot, "current GDP value must be positive"))
	s0 = c0.V ./ tot0
	s1 = c.V ./ tot
	active = (abs.(c0.V) .+ abs.(c.V)) .> 0
	dlnP = 0.0
	for j in eachindex(c.V)
		active[j] || continue
		(c0.V[j] == 0 || c.V[j] == 0) && continue
		(c0.Q[j] <= 0 || c.Q[j] <= 0) && continue
		u0 = abs(c0.V[j]) / c0.Q[j]
		u1 = abs(c.V[j]) / c.Q[j]
		(u0 <= 0 || u1 <= 0) && continue
		dlnP += 0.5 * (s0[j] + s1[j]) * log(u1 / u0)
	end
	return exp(dlnP)
end

"""
	gdp_income(sol::Solution, base::Solution) -> Float64

Income-side real GDP index (ADR-0018): nominal wage-bill growth deflated by
`gdp_deflator`. 1.0 at `base`.
"""
function gdp_income(sol::Solution, base::Solution)::Float64
	c = gdp_components(sol.model, sol)
	c0 = gdp_components(base.model, base)
	return (c.wage_bill / c0.wage_bill) / gdp_deflator(sol, base)
end

"""
	gdp_expenditure(sol::Solution, base::Solution) -> Float64

Expenditure-side Divisia dual (ADR-0018 + ADR-0019): nominal `ΣV` growth
deflated by `gdp_deflator`. No external-financing adjustment is needed: the
transfer F is already booked inside `ΣV` (it enters household expenditure
and hence gross consumption), and the exact identity `wedge = −canary.diff
≈ 0` at η = 1 solutions gives `gdp_income ≡ gdp_expenditure` there —
including fixed F3 cells and the baseline. At BF η = 0 cells the two sides
differ by the factor-market wedge. 1.0 at `base`.
"""
function gdp_expenditure(sol::Solution, base::Solution)::Float64
	c = gdp_components(sol.model, sol)
	c0 = gdp_components(base.model, base)
	return (sum(c.V) / sum(c0.V)) / gdp_deflator(sol, base)
end

"""
	gdp_wedge(sol::Solution) -> Float64

External wedge `ΣV − w·ΣL` at a solved equilibrium (ADR-0018 + ADR-0019):
the negated canary identity gap (`wedge = −canary.diff` exactly); ≈ 0 at
η = 1 solutions, the factor-market gap at BF η = 0 cells. Diagnostic only,
never gated to zero.
"""
function gdp_wedge(sol::Solution)::Float64
	return gdp_components(sol.model, sol).wedge
end

"""
	real_consumption(sol::Solution) -> Float64

Canonical name of the household-consumption (welfare) Törnqvist index stored
in `Solution.real_gdp` (ADR-0018). `real_gdp` itself is unchanged (D3).
"""
function real_consumption(sol::Solution)::Float64
	return sol.real_gdp
end

# --- src/impulses.jl (verbatim) ---

function load_impulses(filename)
  filedir = joinpath(pwd(), "data/", filename)
  impulses = CSV.read(filedir, DataFrames.DataFrame, delim=",", decimal='.', missingstring=["-", "x"]) #read in from csv
  select!(impulses, Not(names(impulses)[1])) # remove the first column
  impulses ./ 1_000_000
end

# --- src/variance_decomposition.jl (verbatim) ---

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
#   - 1 − sum(S_f) is the interaction share for a complete balanced factorial;
#     it need not be zero.  With failed points it also contains missing-data bias.
#
# Author: calculato (AI research assistant)
# ═══════════════════════════════════════════════════════════════════════════════

using ProgressMeter
using Printf
using Statistics
using CSV
using DataFrames

"""
    eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=nothing, verbose=true)

Run the mobile-labor CES model across a vector of η values.
Returns a vector of `Solution` objects.

# Example
```julia
η_grid = [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0]
sols = eta_sweep(data, shocks, 0.5, 0.5, 0.9, η_grid)
```
"""
function eta_sweep(data::Data, shocks::Shocks, θ::Float64, ϵ::Float64, σ::Float64, η_values::AbstractVector{<:Real}; labor_bar::Union{Float64, Nothing}=nothing, verbose::Bool=true)
    η_values = _validate_grids(η_values; eta=true)[1]
    for (name, value) in (("θ", θ), ("ϵ", ϵ), ("σ", σ))
        isfinite(value) || throw(ArgumentError("$name must be finite"))
    end
    N = length(data.factor_share)
    results = Vector{Solution}(undef, length(η_values))

    baseline_init = [ones(N); data.λ; 1.0]
    progress = verbose ? Progress(length(η_values); desc="η sweep: ") : nothing
    try
        for (idx, η) in enumerate(η_values)
            local sol
            try
                model = mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=labor_bar)
                sol = solve(model; init=copy(baseline_init))
            catch e
                throw(ArgumentError("η sweep failed at η=$η, θ=$θ, ϵ=$ϵ, σ=$σ: $(sprint(showerror, e))"))
            end
            results[idx] = sol
            verbose && next!(progress)
        end
    finally
        verbose && finish!(progress)
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
- `values`: the output values at each grid point (non-finite where solver failed)
- `n_failed`: number of grid points with non-finite output (solver failure)
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
On a complete balanced factorial, this uses the standard Sobol/ANOVA
definitions: `S_f = SS_f / SS_total` and `ST_f = 1 - SS_{-f}/SS_total`.
`1 - sum(S_f)` is the interaction share (and is not generally zero).
Failed grid points are retained as NaN and the resulting estimates are
weighted over the valid observations; they should therefore not be treated
as estimates from a balanced design.
`verbose` controls the progress and informational messages.
"""
function variance_decomposition(
    data::Data,
    shocks::Shocks;
    η_values::Vector{Float64} = [0.0, 1.0],
    ϵ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    θ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    σ_values::Vector{Float64} = [0.1, 0.5, 0.99],
    output::Symbol = :real_gdp,
    labor_bar::Union{Float64, Nothing} = nothing,
    verbose::Bool = true,
)

    output in (:real_gdp, :nominal_gdp, :sectoral_q, :sectoral_p) ||
        throw(ArgumentError("unsupported output=$output; choose :real_gdp, :nominal_gdp, :sectoral_q, or :sectoral_p"))
    η_values = _validate_grids(η_values; eta=true)[1]
    ϵ_values, θ_values, σ_values = _validate_grids(ϵ_values, θ_values, σ_values)

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

    baseline_init = [ones(N); data.λ; 1.0]
    progress = verbose ? Progress(n; desc="Variance decomposition: ") : nothing
    try
        for i in 1:n
            η = grid.η[i]
            ϵ = grid.ϵ[i]
            θ = grid.θ[i]
            σ = grid.σ[i]

            try
                model = mobile_labor_model(data, shocks, θ, ϵ, σ, η; labor_bar=labor_bar)
                sol = solve(model; init=copy(baseline_init))

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
            verbose && next!(progress)
        end
    finally
        verbose && finish!(progress)
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

function _validate_grids(grids::AbstractVector{<:Real}...; eta=false)
    checked = Vector{Vector{Float64}}(undef, length(grids))
    for (i, grid) in enumerate(grids)
        isempty(grid) && throw(ArgumentError("factor grid $i must be nonempty"))
        values = Float64.(grid)
        all(isfinite, values) || throw(ArgumentError("factor grid $i must contain only finite values"))
        if eta
            all(v -> v == 0.0 || v == 1.0, values) || throw(ArgumentError(
                "η values must be 0 (immobile) or 1 (fully mobile); intermediate reallocation was retired (ADR-0010)"))
        end
        checked[i] = values
    end
    checked
end


"""
    _compute_sobol_indices(factor_names, grid, values, output_name)

Compute first-order (S_f) and total-order (ST_f) Sobol indices on a
balanced full factorial grid using ANOVA sum-of-squares decomposition.

First-order:   S_f = SS_f / SS_total
Total-order:   ST_f = 1 − SS_{-f} / SS_total

where SS_{-f} is the sum of squares explained by all factors EXCEPT f.

All non-finite values (from solver failures or invalid model output) are
dropped; `n_failed` records how many points were lost. The total-order index
is biased if grid points are missing, so `n_failed` should be reported.
"""
function _compute_sobol_indices(factor_names::Vector{String}, grid::DataFrame, output_vals::Vector{Float64}, output_name::String)
    n_total = length(output_vals)
    valid = isfinite.(output_vals)
    n_valid = sum(valid)
    n_failed = n_total - n_valid
    n_failed > 0 && @warn "Sobol decomposition for $output_name has $n_failed invalid/missing grid values; estimates are not from a complete balanced design"

    if n_valid == 0
        S_f = Dict(f => NaN for f in factor_names)
        ST_f = Dict(f => NaN for f in factor_names)
        return SobolResult(factor_names, S_f, ST_f, grid, output_vals, n_failed, output_name, 0.0, 0.0, NaN)
    end

    v = output_vals[valid]
    grand_mean = mean(v)
    ss_total = sum((v .- grand_mean) .^ 2)

    # Scale-aware tolerance avoids dividing by numerical noise while
    # preserving variation when outputs are merely small in magnitude.
    variance_scale = max(1.0, maximum(abs.(v))^2) * max(1, n_valid)
    zero_tol = 100 * eps(Float64) * variance_scale
    if ss_total <= zero_tol
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
    # For each combination of levels of all other factors, take the mean
    # of y across the valid observations and weight by its actual count.
    ST_f = Dict{String, Float64}()
    for f in factor_names
        other = setdiff(factor_names, [f])
        other_cols = [grid[!, f_name][valid] for f_name in other]
        combos = unique([Tuple([c[i] for c in other_cols]) for i in 1:n_valid])
        ss_except_f = 0.0
        for combo in combos
            mask = trues(n_valid)
            for (j, f_name) in enumerate(other)
                mask = mask .& (other_cols[j] .== combo[j])
            end
            n_in_mask = sum(mask)
            if n_in_mask > 0
                group_mean = mean(v[mask])
                # In an incomplete design use the actual number of valid
                # observations, not the nominal level count.
                ss_except_f += n_in_mask * (group_mean - grand_mean)^2
            end
        end
        ST_f[f] = 1.0 - ss_except_f / ss_total
    end

    # Do not clip substantive out-of-range values: they diagnose an
    # incomplete/unbalanced design. Only remove insignificant roundoff.
    for f in factor_names
        S_f[f] = _normalize_sobol_index(S_f[f], "S_f[$f]")
        ST_f[f] = _normalize_sobol_index(ST_f[f], "ST_f[$f]")
    end

    ss_explained = sum(Base.values(S_f)) * ss_total
    frac_unexplained = 1.0 - sum(Base.values(S_f))

    return SobolResult(factor_names, S_f, ST_f, grid, output_vals, n_failed, output_name, ss_total, ss_explained, frac_unexplained)
end

function _normalize_sobol_index(x::Float64, label::String)
    isfinite(x) || return x
    tol = 100 * eps(Float64) * max(1.0, abs(x))
    abs(x) <= tol && return 0.0
    abs(x - 1.0) <= tol && return 1.0
    (x < 0.0 || x > 1.0) && @warn "$label=$x is outside [0, 1]; retaining diagnostic value"
    x
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
    n_grid = nrow(vd.grid)
    failed_pct = n_grid == 0 ? NaN : 100 * vd.n_failed / n_grid
    println("Grid: $n_grid points, $(vd.n_failed) failed (non-finite) = $(round(failed_pct, digits=1))%")
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
    finite_factors = [f for f in vd.factors if haskey(vd.S_f, f) && isfinite(vd.S_f[f])]
    if isempty(finite_factors)
        println("Dominant factor: unavailable (all first-order indices are non-finite)")
    else
        dom = argmax([vd.S_f[f] for f in finite_factors])
        @printf("Dominant factor: %s (S_f=%.4f)\n", finite_factors[dom], vd.S_f[finite_factors[dom]])
    end

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
    # Only the BF endpoints are kept (ADR-0010); the sweep compares immobile
    # (η = 0) against fully mobile (η = 1).
    η_values = [0.0, 1.0]
    sols = eta_sweep(data, shocks, θ, ϵ, σ, η_values; labor_bar=labor_bar)
    return EtaSweepResult(η_values, sols)
end

"""Run the sweep and return descriptive diagnostics, without a binary decision."""
function eta_sweep_diagnostics(data::Data, shocks::Shocks; θ::Float64=0.5, ϵ::Float64=0.5, σ::Float64=0.9, labor_bar::Union{Float64, Nothing}=nothing)
    esr = eta_sweep_full(data, shocks; θ=θ, ϵ=ϵ, σ=σ, labor_bar=labor_bar)
    variation = vec(std(sectoral_quantities(esr), dims=2))
    vd = variance_decomposition(data, shocks; η_values=[0.0, 1.0],
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
        η_values = [0.0, 1.0],
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
