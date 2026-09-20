
# Plotting entry points.
#
# The plotting implementation lives in the GLMakie package extension
# (`ext/BeyondHultenGLMakieExt.jl`), which Julia loads automatically when
# `GLMakie` and `BeyondHulten` are both loaded in the same session. The
# stubs below keep the exported plotting names defined for headless use and
# raise an instructive error telling the user how to enable plotting.

function _require_glmakie(fname::Symbol)
	error("`BeyondHulten.$(fname)` requires the GLMakie extension, which is not loaded. " *
		"Run `using GLMakie` together with `using BeyondHulten` " *
		"(install it first with `import Pkg; Pkg.add(\"GLMakie\")` if necessary).")
end

"""
	axis_change_in_level!(fig, data, impulses; options)

Left panel of `panel`: bar chart of demand shocks with error bars.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function axis_change_in_level!(args...; kwargs...)
	_require_glmakie(:axis_change_in_level!)
end

"""
	axis_change_in_price!(fig, data, impulse; options)

Right panel of `panel`: price/quantity scatter.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function axis_change_in_price!(args...; kwargs...)
	_require_glmakie(:axis_change_in_price!)
end

"""
	get_color(data, shocks)

Highlight vector for the shocked sectors in `axis_change_in_price!`.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function get_color(args...; kwargs...)
	_require_glmakie(:get_color)
end

"""
	panel(data, impulses; options, name = "panel")

Two-panel impulse-response figure (bar chart + price/quantity scatter).
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function panel(args...; kwargs...)
	_require_glmakie(:panel)
end

"""
	diff_lambda(data, impulses; options, name = "diff_lambda_imp")

Stacked CGE vs. Leontief sectoral decomposition figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function diff_lambda(args...; kwargs...)
	_require_glmakie(:diff_lambda)
end

"""
	effect_of_different_elasticities(shocks, data, gdp_effect_simple; labor_slack_function, name)

Real-GDP elasticity-gradient figure for several elasticity levels.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function effect_of_different_elasticities(args...; kwargs...)
	_require_glmakie(:effect_of_different_elasticities)
end

"""
	comparison_between_labor_slacks(data, shocks, gdp_effect_simple, title)

GDP figure comparing labour-slack specifications.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function comparison_between_labor_slacks(args...; kwargs...)
	_require_glmakie(:comparison_between_labor_slacks)
end

"""
	labor_slack_gradient(data, impulse)

Real-GDP figure along the labour-slack interpolation.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function labor_slack_gradient(args...; kwargs...)
	_require_glmakie(:labor_slack_gradient)
end

"""
	plot_real_gdp_gradient(results; title, cd, leontief, initial, ylims)

Four-panel real-GDP elasticity-gradient figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_real_gdp_gradient(args...; kwargs...)
	_require_glmakie(:plot_real_gdp_gradient)
end

# ── 5x3 evaluation-matrix data layer (headless; WP2 adds the GLMakie methods) ──
#
# Per-cell sectoral data for the executed labour-closure × financing-closure
# matrix (`runs/matrix_5x3-v10-*`): relative (level/baseline) sectoral vectors
# plus the design baseline, tidy frames for plotting, and pure validation
# helpers. Everything here runs without GLMakie; the `plot_matrix_*` entry
# points below are headless stubs that the GLMakie extension implements.

"""
	matrix_labour_order()::Vector{String}

Canonical labour rows of the 5x3 matrix.
"""
function matrix_labour_order()::Vector{String}
	return ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"]
end

"""
	matrix_financing_order()::Vector{String}

Canonical financing columns of the 5x3 matrix.
"""
function matrix_financing_order()::Vector{String}
	return ["F1", "F2", "F3"]
end

"""
	MatrixBaseline(prices, quantities, wages, consumption, imports, exports)

Design baseline levels: raw prices (≈ 1), gross outputs (= `data.λ`),
CPI-deflated real wages (= 1), gross household demand quantities, import
content of final + intermediate demand (value), exports (value).
"""
struct MatrixBaseline
	prices::Vector{Float64}
	quantities::Vector{Float64}
	wages::Vector{Float64}
	consumption::Vector{Float64}
	imports::Vector{Float64}
	exports::Vector{Float64}
end

"""
	MatrixCellData(...)

One matrix cell: identifying strings, sectoral labels, and relative
(level/baseline, 1.0 = unchanged) sectoral vectors, plus the recorded
metrics/diagnostics dicts. Sectors with a zero baseline map to NaN unless
the level is zero too (then 1.0).
"""
struct MatrixCellData
	run_id::String
	labour::String
	financing::String
	eta::Float64
	eta_s::Float64
	status::String
	sectors::Vector{String}
	labels::Vector{String}
	prices::Vector{Float64}
	quantities::Vector{Float64}
	wages::Vector{Float64}
	consumption::Vector{Float64}
	imports::Vector{Float64}
	exports::Vector{Float64}
	metrics::Dict{String,Any}
	diagnostics::Dict{String,Any}
end

"""
	MatrixDataset(design, baseline, cells, provenance, validation)

A full matrix: design name, shared baseline, per-cell data, free-form
provenance, and a per-run validation map (`run_id => "ok"` or a failure
description).
"""
struct MatrixDataset
	design::String
	baseline::MatrixBaseline
	cells::Vector{MatrixCellData}
	provenance::Dict{String,Any}
	validation::Dict{String,String}
end

"""
	_matrix_rel(level, base) -> Vector{Float64}

Elementwise level/baseline ratios: `base_i > 0 ? level_i / base_i :
(level_i == 0 ? 1.0 : NaN)`. Zero-baseline sectors become NaN unless the
level is zero too.
"""
function _matrix_rel(level::AbstractVector, base::AbstractVector)::Vector{Float64}
	n = length(level)
	length(base) == n || throw(DimensionMismatch(
		"level (length $n) and baseline (length $(length(base))) must match"))
	return Float64[base[i] > 0 ? level[i] / base[i] : (level[i] == 0 ? 1.0 : NaN)
		for i in 1:n]
end

"""Expand a scalar wage to the sectoral vector; pass through a sectoral one."""
function _matrix_wages(w::AbstractVector, n::Int)::Vector{Float64}
	length(w) == n && return Float64.(w)
	length(w) == 1 && return fill(Float64(w[1]), n)
	throw(DimensionMismatch("wage vector has length $(length(w)); expected 1 or $n"))
end

"""
	sectoral_trade_flows(model, sol) -> NamedTuple

Per-sector valued trade flows at the equilibrium: `imports_final`
(household + government + investment + programme import content),
`imports_intermediate` (the ADR-0012 intermediate-import leak with the
ADR-0016 CES bill factor), `exports` (domestic sales abroad, no margin).
`imports = imports_final + imports_intermediate`. Mirrors
`external_balance_canary` exactly (same demand block, same bill factor).
"""
function sectoral_trade_flows(model::Model{MobileLaborCES}, sol::Solution)
	p = Float64.(sol.prices_raw)
	y = Float64.(sol.quantities)
	sect = model.options.elasticities.η == 0.0 ||
		model.options.elasticities.eta_s_vec !== nothing
	w = sect ? sol.wages_raw : sol.wages_raw[1]
	blocks = _mobile_market_demand(model, p, y, w; external_transfer = sol.external_transfer)
	m = Float64.(model.data.import_margin)
	c_gross = blocks.c_dom ./ max.(1 .- m, eps(Float64))
	imports_final = p .* m .* (c_gross .+ blocks.additive .+ model.data.gov_demand .+ model.data.exo_demand)
	θ = model.options.elasticities.θ
	ϵ = model.options.elasticities.ϵ
	k = p .^ ϵ .* model.shocks.supply_shock .^ (ϵ - 1) .*
		_intermediate_price(model.data.Ω_raw, p, θ) .^ (1 - ϵ)
	imports_intermediate = k .* (model.data.M_int ./ model.data.λ) .* y
	exports = p .* model.data.exports_demand
	imports = imports_final .+ imports_intermediate
	return (; imports_final = Float64.(imports_final),
		imports_intermediate = Float64.(imports_intermediate),
		exports = Float64.(exports), imports = Float64.(imports))
end

"""
	matrix_baseline(data, ref_sol) -> MatrixBaseline

Design baseline levels: prices/quantities/wages/consumption from the
reference solution; imports/exports from `sectoral_trade_flows`.
"""
function matrix_baseline(data::Data, ref_sol::Solution)::MatrixBaseline
	n = length(data.factor_share)
	flows = sectoral_trade_flows(ref_sol.model, ref_sol)
	return MatrixBaseline(
		Float64.(ref_sol.prices_raw),
		Float64.(ref_sol.quantities),
		_matrix_wages(ref_sol.wages, n),
		Float64.(ref_sol.consumption),
		Float64.(flows.imports),
		Float64.(flows.exports))
end

"""
	matrix_cell_data(run_id, labour, financing, eta, eta_s, status, data, sol, ref_sol;
	    metrics, diagnostics, labels) -> MatrixCellData

Relative (level/baseline) sectoral vectors for one cell. Levels come from
`sol`, the baseline from `matrix_baseline(data, ref_sol)`. `sectors` are the
long IO names; `labels` default to `"1".."N"`.
"""
function matrix_cell_data(run_id::AbstractString, labour::AbstractString,
		financing::AbstractString, eta::Real, eta_s::Real, status::AbstractString,
		data::Data, sol::Solution, ref_sol::Solution;
		metrics::AbstractDict = Dict{String,Any}(),
		diagnostics::AbstractDict = Dict{String,Any}(),
		labels::Union{AbstractVector, Nothing} = nothing)::MatrixCellData
	n = length(data.factor_share)
	base = matrix_baseline(data, ref_sol)
	flows = sectoral_trade_flows(sol.model, sol)
	sectors = String[string(x) for x in data.io.Sektoren[1:n]]
	lab = labels === nothing ? [string(i) for i in 1:n] : String[string(x) for x in labels]
	length(lab) == n || throw(DimensionMismatch(
		"labels has length $(length(lab)); expected $n"))
	return MatrixCellData(
		String(run_id), String(labour), String(financing),
		Float64(eta), Float64(eta_s), String(status),
		sectors, lab,
		_matrix_rel(sol.prices_raw, base.prices),
		_matrix_rel(sol.quantities, base.quantities),
		_matrix_rel(_matrix_wages(sol.wages, n), base.wages),
		_matrix_rel(sol.consumption, base.consumption),
		_matrix_rel(flows.imports, base.imports),
		_matrix_rel(flows.exports, base.exports),
		Dict{String,Any}(metrics), Dict{String,Any}(diagnostics))
end

"""
	matrix_dataset(design, baseline, cells; provenance, validation) -> MatrixDataset
"""
function matrix_dataset(design::AbstractString, baseline::MatrixBaseline,
		cells::AbstractVector{MatrixCellData};
		provenance::AbstractDict = Dict{String,Any}(),
		validation::AbstractDict = Dict{String,String}())::MatrixDataset
	return MatrixDataset(String(design), baseline, Vector{MatrixCellData}(cells),
		Dict{String,Any}(provenance), Dict{String,String}(validation))
end

"""
	matrix_cell_ids(ids, design) -> Vector{String}

Ids that are exactly `<design>-<labour>-<financing>` for the canonical
labour × financing order, filtered to those present in `ids`, canonical
order. Only a trailing `_vN` version suffix of `design` matches `-vN` in run ids
(`matrix_5x3_v10` → `matrix_5x3-v10-…`; `matrix_5x3` and `survey_v2_matrix` unchanged).
"""
function matrix_cell_ids(ids::AbstractVector{<:AbstractString},
		design::AbstractString)::Vector{String}
	prefix = replace(String(design), r"_v(\d+)$" => s"-v\1")
	have = Set(String.(ids))
	out = String[]
	for lab in matrix_labour_order(), fin in matrix_financing_order()
		id = prefix * "-" * lab * "-" * fin
		id in have && push!(out, id)
	end
	return out
end

"""Canonical-order sort key for a cell (unknown ids sort last)."""
function _matrix_cell_key(labour::AbstractString, financing::AbstractString)
	li = findfirst(==(String(labour)), matrix_labour_order())
	fi = findfirst(==(String(financing)), matrix_financing_order())
	return (li === nothing ? typemax(Int) : li, fi === nothing ? typemax(Int) : fi)
end

"""
	matrix_sectoral_frame(ds) -> DataFrame

Tidy long frame: `run_id`, `labour`, `financing`, `sector`, `label`,
`variable`, `rel`. `variable` in `["price", "quantity", "wage",
"consumption", "imports", "exports"]`; rows ordered by canonical cell
order, then variable, then sector index.
"""
function matrix_sectoral_frame(ds::MatrixDataset)::DataFrame
	vars = ("price", "quantity", "wage", "consumption", "imports", "exports")
	order = sortperm(collect(1:length(ds.cells));
		by = i -> _matrix_cell_key(ds.cells[i].labour, ds.cells[i].financing))
	run_id = String[]
	labour = String[]
	financing = String[]
	sector = String[]
	label = String[]
	variable = String[]
	rel = Float64[]
	for i in order
		cell = ds.cells[i]
		n = length(cell.sectors)
		vecs = (cell.prices, cell.quantities, cell.wages,
			cell.consumption, cell.imports, cell.exports)
		for (v, vec) in zip(vars, vecs), j in 1:n
			push!(run_id, cell.run_id)
			push!(labour, cell.labour)
			push!(financing, cell.financing)
			push!(sector, cell.sectors[j])
			push!(label, cell.labels[j])
			push!(variable, v)
			push!(rel, vec[j])
		end
	end
	return DataFrame(run_id = run_id, labour = labour, financing = financing,
		sector = sector, label = label, variable = variable, rel = rel)
end

"""
	matrix_summary_frame(ds) -> DataFrame

One row per cell: `run_id`, `labour`, `financing`, `eta`, `eta_s`,
`status` plus the metrics columns (first-seen order; missing metrics
=> NaN).
"""
function matrix_summary_frame(ds::MatrixDataset)::DataFrame
	mkeys = String[]
	for cell in ds.cells, k in keys(cell.metrics)
		k in mkeys || push!(mkeys, String(k))
	end
	df = DataFrame(run_id = [c.run_id for c in ds.cells],
		labour = [c.labour for c in ds.cells],
		financing = [c.financing for c in ds.cells],
		eta = [c.eta for c in ds.cells],
		eta_s = [c.eta_s for c in ds.cells],
		status = [c.status for c in ds.cells])
	for k in mkeys
		df[!, k] = [haskey(c.metrics, k) ? c.metrics[k] : NaN for c in ds.cells]
	end
	return df
end

"""Default metric keys compared by `validate_cell`."""
const _MATRIX_VALIDATE_METRICS = ("gdp", "gdp_rel", "gdp_expenditure",
	"gdp_expenditure_rel", "gdp_deflator", "gdp_wedge", "consumption",
	"consumption_rel", "employment", "wage", "nominal_gdp",
	"max_abs_price_dev", "external_transfer", "programme_financing",
	"external_position")

"""
	validate_cell(manifest, re_metrics, re_prices, re_quantities, stored_prices, stored_quantities;
	    rtol, atol, price_atol, quantity_atol, metric_keys) -> NamedTuple

Compare a re-solved cell against its recorded run artifacts. Pure (no
disk). Returns `(pass, failures, deltas, skipped)` with `deltas` keyed
`"price"`, `"quantity"`, `"metric:<key>"` and `skipped` listing the metric
keys that were not compared because they are absent from the manifest
and/or the re-solved metrics (each records NaN in `deltas`). A
non-executed status, a non-passing gate, a metric drift beyond
`atol + rtol*|manifest|`, or a stored price/quantity mismatch (or wrong
length) is a failure entry, never an exception. Skipped metrics are never
failures. `stored_prices === nothing` / `stored_quantities === nothing`
skips that check (records NaN in `deltas`).
"""
function validate_cell(manifest::AbstractDict, re_metrics::AbstractDict,
		re_prices, re_quantities, stored_prices, stored_quantities;
		rtol::Real = 1e-8, atol::Real = 1e-8,
		price_atol::Real = 1e-8, quantity_atol::Real = 1e-8,
		metric_keys = _MATRIX_VALIDATE_METRICS)
	failures = String[]
	deltas = Dict{String,Float64}()
	skipped = String[]
	status = get(manifest, "status", "")
	status == "executed" || push!(failures,
		"status is \"$status\", expected \"executed\"")
	gates = get(manifest, "gates", Dict())
	overall = gates isa AbstractDict ? get(gates, "overall", "fail") : "fail"
	overall == "pass" || push!(failures,
		"gates.overall is \"$overall\", expected \"pass\"")
	m_metrics = get(manifest, "metrics", Dict())
	m_metrics = m_metrics isa AbstractDict ? m_metrics : Dict()
	for k in metric_keys
		ks = String(k)
		in_m = haskey(m_metrics, k) || haskey(m_metrics, ks)
		in_r = haskey(re_metrics, k) || haskey(re_metrics, ks)
		mv = haskey(m_metrics, k) ? m_metrics[k] :
			(haskey(m_metrics, ks) ? m_metrics[ks] : nothing)
		rv = haskey(re_metrics, k) ? re_metrics[k] :
			(haskey(re_metrics, ks) ? re_metrics[ks] : nothing)
		if in_m && in_r
			if mv isa Real && rv isa Real
				tol = Float64(atol) + Float64(rtol) * abs(Float64(mv))
				d = abs(Float64(rv) - Float64(mv))
				deltas["metric:" * ks] = d
				d > tol && push!(failures,
					"metric $ks differs: manifest=$mv re-solved=$rv (tol=$tol)")
			else
				deltas["metric:" * ks] = rv == mv ? 0.0 : NaN
				rv == mv || push!(failures,
					"metric $ks differs: manifest=$mv re-solved=$rv")
			end
		else
			deltas["metric:" * ks] = NaN
			push!(skipped, ks)
		end
	end
	if stored_prices === nothing
		deltas["price"] = NaN
	else
		rp = Float64.(collect(re_prices))
		sp = Float64.(collect(stored_prices))
		if length(sp) != length(rp)
			push!(failures,
				"stored prices length $(length(sp)) != re-solved length $(length(rp))")
			deltas["price"] = NaN
		else
			d = isempty(rp) ? 0.0 : maximum(abs, sp .- rp)
			deltas["price"] = d
			d > Float64(price_atol) && push!(failures,
				"stored prices differ: max|stored - re-solved|=$d exceeds atol=$price_atol")
		end
	end
	if stored_quantities === nothing
		deltas["quantity"] = NaN
	else
		rq = Float64.(collect(re_quantities))
		sq = Float64.(collect(stored_quantities))
		if length(sq) != length(rq)
			push!(failures,
				"stored quantities length $(length(sq)) != re-solved length $(length(rq))")
			deltas["quantity"] = NaN
		else
			d = isempty(rq) ? 0.0 : maximum(abs, sq .- rq)
			deltas["quantity"] = d
			d > Float64(quantity_atol) && push!(failures,
				"stored quantities differ: max|stored - re-solved|=$d exceeds atol=$quantity_atol")
		end
	end
	return (; pass = isempty(failures), failures = failures, deltas = deltas, skipped = skipped)
end

"""
	plot_matrix_overview(args...; kwargs...)

5x3 matrix overview figure (one panel per cell).
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_overview(args...; kwargs...)
	_require_glmakie(:plot_matrix_overview)
end

"""
	plot_matrix_wages(args...; kwargs...)

Sectoral-wage panel of the 5x3 matrix figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_wages(args...; kwargs...)
	_require_glmakie(:plot_matrix_wages)
end

"""
	plot_matrix_prices(args...; kwargs...)

Sectoral-price panel of the 5x3 matrix figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_prices(args...; kwargs...)
	_require_glmakie(:plot_matrix_prices)
end

"""
	plot_matrix_quantities(args...; kwargs...)

Sectoral-quantity panel of the 5x3 matrix figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_quantities(args...; kwargs...)
	_require_glmakie(:plot_matrix_quantities)
end

"""
	plot_matrix_consumption(args...; kwargs...)

Sectoral-consumption panel of the 5x3 matrix figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_consumption(args...; kwargs...)
	_require_glmakie(:plot_matrix_consumption)
end

"""
	plot_matrix_trade(args...; kwargs...)

Sectoral-trade (imports/exports) panel of the 5x3 matrix figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_matrix_trade(args...; kwargs...)
	_require_glmakie(:plot_matrix_trade)
end

"""
	save_matrix_figures(args...; kwargs...)

Write the 5x3 matrix figures to disk.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function save_matrix_figures(args...; kwargs...)
	_require_glmakie(:save_matrix_figures)
end
