# ═══════════════════════════════════════════════════════════════════════════════
# src/core/calibration.jl — v3 open-economy calibration (Phase 3)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Promoted from the frozen `cbase2/src/calibration.jl` (Notebook 03b); read
# cbase2 freely, never edit it. The design notes below are preserved from the
# cbase2 header. Status lives only in `registry/closures.toml` (ADR-0003);
# decisions in `docs/decisions/ADR-0006-experiment-entry-point-and-preregistration.md`
# and `docs/decisions/ADR-0010-port-review-fixes.md`.
#
# National-accounts-consistent demand structure (fixes the v2 inconsistency):
#
#   Y = C + I + G + X - M      (expenditure side, domestic output valuation)
#   Household:Disposable income splits into consumption (1-s)E and saving sE.
#   Exogenous injections: government G (tax-financed, F2 balanced),
#   investment I and exports X (backed by saving: S = I + X - M at baseline).
#   Imports M: fixed sector margins m_i on consumption, government,
#   investment and the programme bundle. Exports carry NO margin (domestic
#   sales abroad).
#
# The saving rate is NOT assumed -- it is the data-implied residual
#   s = 1 - sum(c0_gross)/E_h0,  E_h0 = 1 - sum(gG),
# with c0 the baseline household block calibrated so the model clears
# exactly at y = lambda. The identity S = I + X - M then holds at baseline
# up to the clamping/accounting residual (validated in the notebooks).
# The open-economy blocks (gG, inv, expo and the domestic/import split of
# final demand) are read DIRECTLY from the raw IO table stored in `data.io`;
# no processed artifact is consumed and no sector-index argument is needed.
#
# POST-KEYNESIAN ALTERNATIVE (major option, see cbase2/process_comments.md):
# the single s lumps household and corporate saving (Kaldor-style split is
# the natural refinement: s_h on wage income, s_c on profits, two
# households).
#
# ── Sector-exclusion pipeline (2026-09-17 repair, review findings 2 + 3) ──
# `retained_dataset` REPLACES the old `drop_sectors` matrix slicer. Rather
# than slicing the already-derived full-table matrices (which left Ω rows
# summing to 0.971 (70s) / 0.549 (reduced) and λ, labor_share and the GDP
# aggregates in mixed units), it removes the dropped sectors' rows and
# columns from the RAW IO table and re-runs the full §4.1 accounting
# transformation on the slice:
#
#   * each retained user's Ω_raw row is renormalized on the retained
#     suppliers, i.e. the dropped suppliers' input share is redistributed
#     proportionally over the kept ones (modelling assumption: retained
#     inputs are proportional substitutes for the dropped ones);
#   * λ, factor_share, labor_share and gross output are re-derived on the
#     retained economy, so Σ labor_share = 1 and E_h0 = 1 - ΣgG is the
#     correct baseline-income normalization;
#   * the dropped sectors' final demand leaves the accounting with their
#     rows; their proportional import/product-tax content in the (category-
#     level) FD import and tax rows is removed by the same share, so the
#     domestic/import split of the retained categories is preserved.
#
# The construction asserts Σ_j Ω_raw[u, j] = 1 and Σ labor_share = 1, so the
# CES/CD price-index contract and the income-unit normalization hold for every
# dataset variant instead of silently contaminating the θ-continuation.
# ---------------------------------------------------------------------------

# ------------------------------------------------------------------
# Dataset variants (user-directed design, Notebook 03b).
#   "full"    all 71 sectors (reference; carries the self-loop instability)
#   "70s"     sector 71 dropped (the documented decision; current working set)
#   "reduced" additionally drops every sector with self-share > 0.45
#             (48, 18, 19, 53, 58, 13, 68) -- the genuine-data class proven
#             by the (f) diagnostic; a modelling choice, NOT a data
#             correction. The instability-vs-coverage trade-off is measured
#             by comparing results across variants (notebook 03b).
# Variants are rebuilt from the raw IO table by `retained_dataset`; the
# coverage cost is reported by `dataset_coverage` (pass the FULL-table Data).
# ------------------------------------------------------------------
const DATASET_VARIANTS = Dict{String,Vector{Int}}(
	"full"    => Int[],
	"70s"     => [71],
	"reduced" => [71, 48, 18, 19, 53, 58, 13, 68],
)

"""
	retained_io_table(io, drops; number_sectors = 71) -> DataFrame

Slice the raw IO table down to the retained sectors: dropped sector rows and
columns are removed, the remaining intermediate/VA cells are untouched. The
label column and all non-sector columns (final demand, imports, taxes, value
added, aggregates) are kept.

Imports and goods taxes on final demand are reported BY CATEGORY, not by
supplying sector (the §4.1 source limitation). The dropped sectors' final
demand therefore cannot be cut out of those rows cell-by-cell, so each FD
category's import and goods-tax entries are scaled by the retained share of
that category's purchaser-price final demand. This removes the dropped
products' proportional import/tax content and preserves the domestic/import
split of the retained categories instead of reallocating all category imports
onto fewer products.

The resulting table has the retained sector rows first (in ascending order),
so `generate_data(io; number_sectors = length(keep))` can process it like a
standalone table. The source table's aggregate rows/columns (totals, valuation
rows) are carried along verbatim and are not used by `generate_data`; the
authoritative retained-economy fields are those recomputed by
`retained_dataset`.
"""
function retained_io_table(io::DataFrames.DataFrame, drops::AbstractVector{<:Integer};
		number_sectors::Integer = 71)
	number_sectors >= 1 || throw(ArgumentError("number_sectors must be positive"))
	keep = setdiff(1:Int(number_sectors), Int.(drops))
	isempty(keep) && throw(ArgumentError("cannot drop all sectors"))
	all(d -> 1 <= d <= number_sectors, drops) ||
		throw(ArgumentError("drop indices must lie in 1:number_sectors"))
	# Sector rows are 1:number_sectors, sector columns follow the label column.
	colsel = vcat(1, (1 .+ keep), collect((number_sectors + 2):size(io, 2)))
	rowsel = vcat(keep, collect((number_sectors + 1):size(io, 1)))
	reduced = io[rowsel, colsel]

	# ── Proportional removal of the dropped sectors' FD imports/product taxes ──
	FD_full = final_demand_columns(io)
	FD_red = final_demand_columns(reduced)
	r_imp_full = findfirst(==("Verwendung der Importe"), io.Sektoren)
	r_tx_full  = findfirst(==("Gütersteuern abzüglich Gütersubventionen"), io.Sektoren)
	r_imp_red = findfirst(==("Verwendung der Importe"), reduced.Sektoren)
	r_tx_red  = findfirst(==("Gütersteuern abzüglich Gütersubventionen"), reduced.Sektoren)
	(r_imp_full === nothing || r_tx_full === nothing ||
	 r_imp_red === nothing || r_tx_red === nothing) && throw(ArgumentError(
		"IO table must contain the \"Verwendung der Importe\" and " *
		"\"Gütersteuern abzüglich Gütersubventionen\" rows"))
	n_keep = length(keep)
	for j in FD_red
		# The source table stores integers; the proportional scaling below needs
		# float cells (the FD columns are read as Float64 again by generate_data).
		reduced[!, j] = Float64.(reduced[!, j])
	end
	for k in eachindex(FD_full)
		fd_full = sum(io[1:number_sectors, FD_full[k]])
		fd_keep = sum(reduced[1:n_keep, FD_red[k]])
		share = fd_full > 0 ? fd_keep / fd_full : 0.0
		reduced[r_imp_red, FD_red[k]] *= share
		reduced[r_tx_red, FD_red[k]] *= share
	end
	return reduced
end

"""
	retained_dataset(data_full, drops) -> Data

Rebuild the v1 accounting-consistent dataset for the retained sectors directly
from the raw IO table:

1. slice the raw table to the retained sector rows/columns
   (`retained_io_table`);
2. re-run `generate_data` on the slice, so `Ω_raw` rows are normalized on the
   retained suppliers and `factor_share`, `λ`, `labor_share`, gross output and
   the GDP aggregates are recomputed on the retained economy (the dropped
   sectors' value added and final demand leave with their rows);
3. assert the normalization contracts loudly: every `Ω_raw` row sums to 1 and
   Σ `labor_share` = 1 (baseline factor income at w = 1 is the model's income
   unit).

The v3 open-economy fields are left at their closed-absorption defaults; call
`recalibrate_open(data; exo_scale = ...)` afterwards. `drops` are indices into
`data_full`'s sector list; `drops == []` returns `data_full` unchanged (the
full-table baseline is not rebuilt). `household_baseline` keeps the legacy
Törnqvist default of `assemble_data` (ADR-0005).
"""
function retained_dataset(data_full::Data, drops::Vector{Int})
	n_full = length(data_full.factor_share)
	all(d -> 1 <= d <= n_full, drops) || throw(ArgumentError(
		"drop indices must lie in 1:n_full (data_full has $n_full sectors); got $drops"))
	keep = setdiff(1:n_full, drops)
	isempty(keep) && throw(ArgumentError("cannot drop all sectors"))
	length(keep) == n_full && return data_full
	io = retained_io_table(data_full.io, drops; number_sectors = n_full)
	n = length(keep)
	d = generate_data(io; number_sectors = n)

	# ── Normalization contracts (fail loudly instead of contaminating the θ path) ──
	rowsum = vec(sum(d.Ω_raw; dims = 2))
	all(x -> isapprox(x, 1.0; rtol = 1e-10, atol = 1e-12), rowsum) ||
		error("retained_dataset: Ω_raw rows must sum to 1 after the rebuild " *
		      "(min = $(minimum(rowsum)), max = $(maximum(rowsum)))")
	isapprox(sum(d.labor_share), 1.0; rtol = 1e-10, atol = 1e-12) ||
		error("retained_dataset: Σ labor_share must be 1 after the rebuild " *
		      "(got $(sum(d.labor_share)))")

	return assemble_data(io, d)
end

"""
	dataset_coverage(data_v1, drops) -> NamedTuple

Coverage report for a sector-drop configuration: shares of gross output,
value added and household final demand kept, plus the dropped sector indices.
Pass the FULL-table `Data` (the shares are measured against the full economy);
used by notebook 03b to document the "reduced" variant's cost.
"""
function dataset_coverage(data_v1::Data, drops::Vector{Int})
	dropped = sort(unique(drops))
	g = sum(data_v1.gross_output_basic)
	v = sum(data_v1.value_added)
	f = sum(data_v1.domestic_final_demand)
	gk = g - sum(data_v1.gross_output_basic[dropped])
	vk = v - sum(data_v1.value_added[dropped])
	fk = f - sum(data_v1.domestic_final_demand[dropped])
	(; dropped_sectors = dropped,
	   gross_share_kept = round(100 * gk / g; digits=2),
	   va_share_kept = round(100 * vk / v; digits=2),
	   fd_share_kept = round(100 * fk / f; digits=2))
end

"""
	recalibrate_open(data; exo_scale = 1.0) -> Data

v3 open-economy recalibration of a v1 `Data` object, derived directly from the
raw IO table stored in `data.io` (the former `cbroot` artifact
AC_domestic_final_demand.csv is not needed: it is reproduced exactly from the
table by `final_demand_split`). No `drops` argument: `retained_dataset` has
already rebuilt the table for the retained sectors, and the open-economy
blocks are read from that table. Returns a new Data with:
  gov_demand      gross government consumption (margin applies)
  exo_demand      gross investment: equipment + construction + inventories
  exports_demand  exports (no margin)
  household_baseline  c0_gross: residual baseline household demand
  consumption_share   omega = c0_gross / sum(c0_gross)  (CPI weights, sum 1)
  import_margin       sector margin m_i of total final demand
  saving_rate         s = 1 - sum(c0_gross)/(1 - sum(gG))
Asserts the finiteness gate: round-gain column sums strictly below 1.
"""
function recalibrate_open(data::Data; exo_scale::Real = 1.0)
	(; Ω_raw, factor_share, λ, io) = data
	n = length(factor_share)
	fs = factor_share
	0 <= exo_scale <= 1 || throw(ArgumentError("exo_scale must be in [0, 1]"))

	# ── Final demand by sector and category, domestic split (raw table) ──
	# Categories (final_demand_split order): 1 private_consumption,
	# 2 private_orgs, 3 government, 4 equipment, 5 construction,
	# 6 inventories, 7 exports.
	fd = final_demand_split(io, n)
	tot = fd.tot
	dom = fd.dom
	scale = 1.0 / data.gdp_production          # EUR m -> model units (GDP_P = 1)

	# Sector import margin of total final demand (proportional §4.1 split)
	m = 1.0 .- vec(sum(dom; dims=2)) ./ max.(vec(sum(tot; dims=2)), 1e-12)
	m = clamp.(m, 0.0, 1.0)

	gG = tot[:, 3] .* scale                        # government (gross, margin applies)
	inv = (tot[:, 4] .+ tot[:, 5] .+ tot[:, 6]) .* scale .* exo_scale  # investment (gross, margin applies)
	expo = tot[:, 7] .* scale .* exo_scale         # exports (NO margin -- domestic sales abroad)

	# ── Baseline household block: the table's own domestic household final ──
	# demand (LaForge exact A-bill fix, 2026-09-17). The intermediate demand is
	# charged with the DOMESTIC bill A_u (row 73) distributed by the domestic
	# technology Ω_raw — NOT the purchaser-price bill (1−fs)·λ ≡ A + imports +
	# taxes. By the table's row identity the residual c0 is then exactly the
	# observed domestic household final demand: ZERO negative sectors, no
	# clamping. The imported (row 74) and taxed (row 75) intermediate content is
	# an explicit external-account leak carried in data.M_int / data.T_int
	# (ADR-0012, ADR-0013); leaving row 75 out breaks the external-account canary
	# by exactly that term.
	Mλ = Ω_raw' * data.A_bill
	c0_dom = λ .- Mλ .- (1.0 .- m) .* (gG .+ inv) .- expo
	# Floating-point dust is clamped; genuine negative MASS fails. Measured:
	# full-71 is exact (dust ~4e-19); 70s carries ONE microscopic negative
	# (−2.4e-6, the retained-economy identity gap — LaForge caution), clamped
	# and documented; "reduced" would fail its identity gap (~6e-3) and is
	# deferred as a variant (dataset-variant decision pending).
	dust = -sum(min.(c0_dom, 0.0))
	dust > 1e-5 && error("A-bill household residual negative mass ", dust,
		" in ", count(<(0), c0_dom), " sectors: the domestic-bill identity is broken")
	c0_dom = max.(c0_dom, 0.0)
	c0_gross = c0_dom ./ max.(1.0 .- m, 1e-6)

	# ── Saving rate (data-implied) and CPI weights ──
	E_h0 = 1.0 - sum(gG)                # baseline after-tax income (baseline income = 1)
	s = 1.0 - sum(c0_gross) / E_h0
	# s < 0 admissible inside the exo range: shrinking injections raise the
	# residual household demand above income — the excess is externally
	# financed (S = I + X − M with negative S), never clamped.
	@assert -1.0 < s < 1.0 "calibrated saving rate out of range: $s"
	ω = c0_gross ./ sum(c0_gross)

	# ── Finiteness gate: worst-case round-gain column sums (F2, marginal tau = 0) ──
	# Direct intermediate round: the DOMESTIC bill coefficient a_u = A_bill/λ_u
	# (replaces the old total-bill coefficient 1−fs_u); plus the consumption
	# round through labour income. Enforced only for s >= 0.
	a = data.A_bill ./ λ
	colsums = a .+ (1.0 .- m) .* (1.0 - s) .* fs
	s >= 0 && @assert all(<(1), colsums) "finiteness gate failed: max column sum = $(maximum(colsums))"

	println("recalibrate_open (v3, N=", n, "): saving rate s = ", round(s; digits=4),
		", tau0 = ", round(sum(gG); digits=4),
		", export share = ", round(sum(expo); digits=4),
		", investment share = ", round(sum(inv); digits=4),
		", intermediate imports = ", round(sum(data.M_int); digits=4))

	return Data(data.io, data.Ω, data.Ω_raw, ω, data.factor_share, data.λ,
		data.labor_share, data.consumption_share_gross_output, data.grossy,
		data.value_added, data.gross_output_basic, data.value_added_components,
		data.imports_intermediate, data.import_share, data.domestic_final_demand,
		gG, c0_gross, m, inv, expo, s,
		data.A_bill, data.M_int, data.T_int,
		data.gdp_production, data.gdp_income, data.gdp_expenditure)
end
