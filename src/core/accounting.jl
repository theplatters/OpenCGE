# src/core/accounting.jl — types, Data, Shocks, Model, generate_data, read_data
# (relocated verbatim from src/interface.jl; closure types live in src/closures/labor/types.jl)

abstract type AbstractElasticities end
abstract type AbstractData end

"""Abstract supertype for financing-closure descriptions (concrete types live in
`src/closures/financing/financing.jl`)."""
abstract type AbstractFinancing end


"""
	CESElasticities(θ::Float64, ϵ::Float64, σ::Float64)
elasticity parameters for the CES model. θ = Elasticity between goods, ϵ = elasticity of substitution between labor and goods, σ = elasticity of consumption
"""
struct CESElasticities <: AbstractElasticities
	θ::Float64
	ϵ::Float64
	σ::Float64
end

struct CobbDouglasElasticities <: AbstractElasticities
	α::Vector{Float64}
	β::Vector{Float64}
end

struct LeontiefElasticies <: AbstractElasticities
end

struct LeontiefElasticiesLabor <: AbstractElasticities
end

abstract type ModelType end
struct CES <: ModelType
	elasticities::CESElasticities
	labor_slack::Union{Function, Symbol}
	labor_reallocation::Bool
end

CES() = CES(CESElasticities(0.01, 0.5, 0.9), full_labor_slack, false)
CES(elasticities::CESElasticities, labor_slack) = CES(elasticities, labor_slack, false)
CES(elasticities::CESElasticities; labor_slack=full_labor_slack, labor_reallocation=false) =
	CES(elasticities, labor_slack, labor_reallocation)
struct Leontief <: ModelType
	labor_effect::Bool
end

Leontief() = Leontief(true)

struct CobbDouglas <: ModelType
	elasticities::CobbDouglasElasticities
	labor_slack::Union{Function, Symbol}
end


struct Data <: AbstractData
	io::DataFrame
	# Conditional intermediate-input shares in B&F orientation (DOMESTIC only):
	# rows are using sectors and columns are supplying sectors. Each row sums to 1.
	Ω::Matrix{Float64}
	# Raw (domestic + imported) conditional input-share matrix, kept for audit.
	Ω_raw::Matrix{Float64}
	consumption_share::Vector{Float64}
	factor_share::Vector{Float64}
	λ::Vector{Float64}                       # standard Domar weights (gross_output / GDP)
	labor_share::Vector{Float64}            # composite VA weight = λ .* factor_share
	consumption_share_gross_output::Vector{Float64}
	grossy::Vector{Float64}                 # gross output at basic prices (= Produktionswert)
	value_added::Vector{Float64}            # gross value added at basic prices
	# --- §4.1 accounting-consistent extensions ---
	gross_output_basic::Vector{Float64}
	value_added_components::DataFrames.DataFrame
	imports_intermediate::Vector{Float64}
	import_share::Vector{Float64}
	domestic_final_demand::Vector{Float64}
	# ── v3 open-absorption calibration (ported from cbase2, Notebook 03b) ──
	# Field names/order match cbase2/src/core/interface.jl exactly. Compatibility
	# defaults (compact constructor, read_data) are zeros — EXCEPT
	# household_baseline, which defaults to the legacy Törnqvist base
	# (consumption_share .* sum(labor_share)) so legacy real_gdp behavior is
	# preserved; the Phase 3 calibration overwrites it (ADR-0005).
	gov_demand::Vector{Float64}        # exogenous government consumption (gross, margin applies)
	household_baseline::Vector{Float64} # baseline household block (Törnqvist base; legacy default = cs .* Σ labor_share)
	import_margin::Vector{Float64}     # sector import share of MARGINAL final demand (0 = closed absorption)
	exo_demand::Vector{Float64}        # exogenous investment: equipment + construction + inventories (margin applies)
	exports_demand::Vector{Float64}    # exogenous exports (domestic sales abroad -- NO import margin)
	saving_rate::Float64               # calibrated non-consumption share of after-tax income (S = I + X - M)
	# ── A-bill calibration (LaForge exact domestic-bill fix, 2026-09-17) ──
	# The model's intermediate demand is charged with the DOMESTIC intermediate
	# bill (raw-table row 73, "Gesamte Verwendung der inländischen Produktion"),
	# NOT the purchaser-price bill (1−fs)·λ ≡ A + imports + taxes (rows 74/75).
	# The imported+taxed content is an explicit external-account leak (M_int).
	A_bill::Vector{Float64}            # per-user domestic intermediate bill at baseline (model units)
	M_int::Vector{Float64}             # per-user intermediate imports at baseline (model units)
	gdp_production::Float64
	gdp_income::Float64
	gdp_expenditure::Float64
end

"""
	conditional_input_shares(intermediate_use)

Convert an input-output block from the Destatis convention to the conditional
intermediate-input share matrix used in the Baqaee--Farhi replication code.

`intermediate_use[s, u]` is the value supplied by sector `s` to using sector
`u`. The returned matrix is oriented as user-by-supplier and has elements
`Ω[u, s] = intermediate_use[s, u] / sum_r(intermediate_use[r, u])`, so every
row sums to one. The outer intermediate share `(1 - factor_share[u])` is
applied separately by the model.
"""
function conditional_input_shares(intermediate_use::AbstractMatrix{<:Real})
	user_by_supplier = permutedims(intermediate_use)
	intermediate_totals = sum(user_by_supplier, dims = 2)
	any(iszero, intermediate_totals) &&
		throw(ArgumentError("every using sector must have positive intermediate use"))
	return user_by_supplier ./ intermediate_totals
end

function Data(filename::String)
	read_data(filename)
end

# Compact constructor retained for small synthetic fixtures and older clients.
# The accounting-extension fields are immaterial for equilibrium calculations.
# Compatibility defaults: v3 absorption fields are zeros, EXCEPT household_baseline
# which defaults to the legacy Törnqvist base (consumption_share .* Σ labor_share),
# NOT cbase2's zeros default — so legacy real_gdp behavior is preserved (ADR-0005).
function Data(io::DataFrame, Ω::AbstractMatrix, consumption_share::AbstractVector,
		factor_share::AbstractVector, λ::AbstractVector, labor_share::AbstractVector,
		consumption_share_gross_output::AbstractVector, grossy::AbstractVector,
		value_added::AbstractVector)
	n = length(grossy)
	all(length(v) == n for v in (consumption_share, factor_share, λ, labor_share,
		consumption_share_gross_output, value_added)) ||
		throw(DimensionMismatch("synthetic Data vectors must have a common length"))
	size(Ω) == (n, n) ||
		throw(DimensionMismatch("synthetic Data Ω must have size (n, n)"))
	Ωf = Matrix{Float64}(Ω)
	gv = Float64.(grossy)
	va = Float64.(value_added)
	cs = Float64.(consumption_share)
	ls = Float64.(labor_share)
	fs_f = Float64.(factor_share)
	λ_f = Float64.(λ)
	household_baseline = cs .* sum(ls)
	# Closed-synthetic default: with no import data the DOMESTIC bill equals
	# the total bill, so a_u = 1−fs_u and the demand block is identical to the
	# pre-A-bill behaviour (the goldens and legacy clients keep their values).
	A_bill = (1.0 .- fs_f) .* λ_f
	Data(io, Ωf, Ωf, cs, fs_f, λ_f,
		ls, Float64.(consumption_share_gross_output), gv, va, gv, DataFrame(),
		zeros(n), zeros(n), zeros(n), zeros(n), household_baseline, zeros(n), zeros(n),
		zeros(n), 0.0, A_bill, zeros(n), sum(va), sum(va), sum(va))
end

# ── A-bill extraction (LaForge exact domestic-bill fix, 2026-09-17) ──
# Row 73 "Gesamte Verwendung der inländischen Produktion" = each using sector's
# DOMESTIC intermediate bill; row 74 "Verwendung der Importe" = its imported
# intermediates. Both scaled to model units by 1/GDP_P. Located by label so the
# rebuilt retained-sector tables (shorter sector block) work identically.
function _domestic_bills(io::DataFrames.DataFrame, n::Int, gdp::Real)
	lbl = String.(io.Sektoren)
	rA = findfirst(x -> occursin("inländischen Produktion", x), lbl)
	rImp = findfirst(==("Verwendung der Importe"), lbl)
	(rA === nothing || rImp === nothing) && throw(ArgumentError(
		"IO table must contain the domestic-use row (73) and the import-use row (74)"))
	A_bill = [Float64(io[rA, 1 + u]) for u in 1:n] ./ gdp
	M_int = [Float64(io[rImp, 1 + u]) for u in 1:n] ./ gdp
	return A_bill, M_int
end

# Final-demand categories of the Destatis table, in table order. They are
# looked up by NAME (not by the full-table position 75:81) so the accounting
# transformation also works on rebuilt retained-sector tables, where the
# sector block is shorter.
const _FD_COLUMN_NAMES = [
	"Konsumausgaben der privaten Haushalte im Inland",
	"Konsumausgaben der privaten Organisationen o.E.",
	"Konsumausgaben des Staates",
	"Anlageinvestitionen f.Ausrüstungen u.sonst.Anlagen",
	"Anlageinvestitionen für Bauten",
	"Vorratsveränderungen und Nettozugang an Wertsachen",
	"Exporte",
]

"""
	final_demand_columns(io) -> Vector{Int}

Column positions of the seven final-demand categories in `io`, located by
name so the function is robust to rebuilt retained-sector tables (shorter
sector block) and to column reordering.
"""
function final_demand_columns(io::DataFrames.DataFrame)
	FD = Int[]
	for name in _FD_COLUMN_NAMES
		j = findfirst(==(name), names(io))
		j === nothing && throw(ArgumentError(
			"IO table is missing the final-demand column \"$name\""))
		push!(FD, j)
	end
	return FD
end

"""
	final_demand_split(io, number_sectors) -> NamedTuple

Final-demand blocks at purchaser prices for the first `number_sectors` sector
rows of `io`, plus the §4.1 proportional domestic/import/tax split:

  tot[s, k]       purchaser-price final demand of sector s, category k
  dom[s, k]       domestic content (the rest is supplied by imports + taxes)
  imp_final[k]    total imported final demand of category k
  imptx_final[k]  goods taxes less subsidies on category k

The domestic share of each category is uniform across sectors (the source
table reports imports and product taxes by category, not by supplying
sector).
"""
function final_demand_split(io::DataFrames.DataFrame, number_sectors::Integer)
	FD = final_demand_columns(io)
	r_imp = findfirst(==("Verwendung der Importe"), io.Sektoren)
	r_tx  = findfirst(==("Gütersteuern abzüglich Gütersubventionen"), io.Sektoren)
	(r_imp === nothing || r_tx === nothing) && throw(ArgumentError(
		"IO table must contain the \"Verwendung der Importe\" and " *
		"\"Gütersteuern abzüglich Gütersubventionen\" rows"))
	tot = Matrix{Float64}(io[1:number_sectors, FD])
	imp_final = vec(Matrix{Float64}(io[r_imp:r_imp, FD]))
	imptx_final = vec(Matrix{Float64}(io[r_tx:r_tx, FD]))
	dom = similar(tot)
	for k in eachindex(FD)
		cat_total = sum(tot[:, k])
		domfrac = cat_total > 0 ? (cat_total - imp_final[k] - imptx_final[k]) / cat_total : 0.0
		dom[:, k] = tot[:, k] .* max(domfrac, 0.0)
	end
	return (; FD, tot, dom, imp_final, imptx_final)
end

"""
	generate_data(io::DataFrames.DataFrame; number_sectors::Int = 71)

Pulls the key econometric variables used by the model out of the extended IO
table and performs the §4.1 accounting-consistency transformation:

1. separates domestic vs imported intermediate and final uses;
2. decomposes value added into wages / other production taxes / depreciation /
   net operating surplus (never relabelling total VA as "labour");
3. computes gross output at basic prices and the composite factor share;
4. applies a guard when the three GDP sides differ by 10% or more;
5. builds the raw conditional input-share matrix `Ω_raw` used by equilibrium
   technology, while retaining domestic `Ω` for audit and incidence data.

`number_sectors` is the number of modeled sector rows/columns in `io` (71 for
the full Destatis table). The retained-sector pipeline (`retained_dataset` in
`src/core/calibration.jl`) rebuilds the IO table without the dropped sectors
and passes the reduced count, so `Ω_raw`, `λ`, `labor_share` and the GDP
aggregates are recomputed on the retained economy: every `Ω_raw` row is
normalized on the retained suppliers and Σ `labor_share` is exactly 1,
instead of inheriting full-table shares (the `drop_sectors` defect).

Returns a NamedTuple; `assemble_data` binds it to the table in a `Data`
object.
"""
function generate_data(io::DataFrames.DataFrame; number_sectors::Int = 71)
	number_sectors >= 1 || throw(ArgumentError("number_sectors must be positive"))
	SEC = 2:number_sectors+1          # sector columns (supplier rows / user cols)

	# --- Imports and product taxes are reported BY COLUMN (using sector + final cat).
	# Index those rows DIRECTLY by their table position (rows are separate from the
	# sector block); the final-demand block is located by name. ---
	r_imp = findfirst(==("Verwendung der Importe"), io.Sektoren)
	r_tx  = findfirst(==("Gütersteuern abzüglich Gütersubventionen"), io.Sektoren)
	(r_imp === nothing || r_tx === nothing) &&
		throw(ArgumentError("IO table is missing the import / goods-tax rows"))
	fd = final_demand_split(io, number_sectors)
	imp_inter_byuser = Matrix{Float64}(io[r_imp:r_imp, SEC])[:]
	imp_final        = fd.imp_final
	imptx_final      = fd.imptx_final

	# Raw intermediate-use matrix Z[s,u] (domestic + imported, combined).
	Z = Matrix{Float64}(io[1:number_sectors, SEC])

	# Imported final demand by category.
	imp_final_total  = sum(imp_final)
	imptx_final_total = sum(imptx_final)

	# Domestic intermediate matrix: allocate imports proportionally across suppliers,
	# preserving each user's domestic supplier composition (stated assumption; the
	# source table reports imports only by using sector, not by supplying sector).
	inter_into_user     = sum(Z; dims=1)[:]
	dom_inter_into_user = inter_into_user .- imp_inter_byuser
	denom = max.(inter_into_user, 1e-12)
	# Guard: the source table records imports exceeding intermediate use for a few
	# sectors (a known per-sector data inconsistency). Clamp the domestic remainder
	# to zero so the domestic intermediate matrix stays non-negative; those sectors
	# are then treated as fully import-dependent in the domestic technology.
	dom_inter_into_user = max.(dom_inter_into_user, 0.0)
	Z_dom = Z .* (dom_inter_into_user' ./ denom')
	@assert isapprox(sum(Z_dom; dims=1)[:], dom_inter_into_user; rtol=1e-9)

	# Conditional DOMESTIC input-share matrix (user-by-supplier); each user row sums to 1.
	# Normalize by the per-USER (column) total: divide Z_dom' (rows = users) by the
	# column vector of user totals so every row sums to one.
	Ω_dom = Z_dom' ./ max.(vec(sum(Z_dom; dims=1)), 1e-12)
	# Raw (domestic + imported) share matrix, kept for audit/reference.
	Ω_raw = conditional_input_shares(Z)

	# --- Value added: decompose; do NOT relabel total VA as "labour" ---
	gva     = Vector{Float64}(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), SEC])
	wage    = Vector{Float64}(io[findfirst(==("Arbeitnehmerentgelt im Inland"), io.Sektoren), SEC])
	othertx = Vector{Float64}(io[findfirst(==("Sonst.Produktionsabgaben abzgl. sonst.Subventionen"), io.Sektoren), SEC])
	dep     = Vector{Float64}(io[findfirst(==("Abschreibungen"), io.Sektoren), SEC])
	netop   = Vector{Float64}(io[findfirst(==("Nettobetriebsüberschuss"), io.Sektoren), SEC])
	prodval = Vector{Float64}(io[findfirst(==("Produktionswert"), io.Sektoren), SEC])
	va_total = wage .+ othertx .+ dep .+ netop
	@assert isapprox(va_total, gva; rtol=1e-9) "VA components must sum to gross value added"

	# --- Gross output / value added at basic prices ---
	grossy = prodval                                  # basic-price gross output
	value_added = gva
	factor_share = gva ./ grossy                      # composite VA share (NOT pure labour)
	import_share = imp_inter_byuser ./ max.(grossy, 1e-12)

	# --- GDP three-side reconciliation (basic prices) ---
	GDP_P = sum(gva)
	GDP_I = sum(wage) + sum(othertx) + sum(dep) + sum(netop)
	fd_purch_total = sum(fd.tot)
	fd_dom_basic   = fd_purch_total - imp_final_total - imptx_final_total
	GDP_E = fd_dom_basic
	@assert GDP_P ≈ GDP_I "production must equal income"
	res_P_E = abs(GDP_P - GDP_E) / GDP_P
	println("GDP reconciliation: P=$(round(GDP_P))  I=$(round(GDP_I))  E=$(round(GDP_E))  ",
	        "fd_purch=$(round(fd_purch_total))  imp_final=$(round(imp_final_total))  ",
	        "imptx_final=$(round(imptx_final_total))  residual=$(round(res_P_E*100; digits=3))%")
	@assert res_P_E < 0.10 "production vs expenditure discrepancy > 10% (gross accounting error)"

	# --- Final-demand absorption (B&F orientation) ---
	# Household absorption is calibrated against the same raw technology matrix
	# used by equilibrium production. Ω_dom remains available for domestic
	# incidence/audit calculations.
	consumption_share = (I - Diagonal(1.0 .- factor_share) * Ω_raw)' * grossy
	@views consumption_share[consumption_share .< 0] .= 0
	consumption_share = consumption_share / sum(consumption_share)

	# Standard Domar weights: λ_s = gross_output_s / GDP (always positive).
	λ = grossy ./ GDP_P

	# Model's historical "labour share" = composite VA weight (explicit, not relabelled).
	labor_share = λ .* factor_share

	# --- Domestic final-demand vector at basic prices (shock incidence target) ---
	# `fd.dom` is the §4.1 proportional domestic split of the purchaser-price
	# categories (see `final_demand_split`).
	domestic_final_demand = vec(sum(fd.dom; dims=2))
	consumption_share_gross_output = domestic_final_demand ./ grossy

	va_comp = DataFrame(
		sector = io.Sektoren[1:number_sectors],
		wage = wage,
		other_production_tax = othertx,
		depreciation = dep,
		net_operating_surplus = netop,
		gross_value_added = gva,
	)

	return (
		Ω = Ω_dom,
		Ω_raw = Ω_raw,
		consumption_share = consumption_share,
		factor_share = factor_share,
		λ = λ,
		labor_share = labor_share,
		consumption_share_gross_output = consumption_share_gross_output,
		grossy = grossy,
		value_added = value_added,
		gross_output_basic = prodval,
		value_added_components = va_comp,
		imports_intermediate = imp_inter_byuser,
		import_share = import_share,
		domestic_final_demand = domestic_final_demand,
		gdp_production = GDP_P,
		gdp_income = GDP_I,
		gdp_expenditure = GDP_E,
	)
end


"""
	Shocks

Store equal-length sectoral supply, household-demand, raw-demand, autonomous-demand,
and investment-demand vectors. Use the keyword constructor when the two additive
final-demand vectors should retain their zero defaults.

# Example
```julia-repl
julia> Shocks(ones(length(data.grossy)), ones(length(data.grossy)))
```
"""
struct Shocks
	supply_shock::Vector{Float64}
	demand_shock::Vector{Float64}
	demand_shock_raw::Vector{Float64}
	# Autonomous (extra-household) final demand, e.g. exports or government
	# spending. The mobile model scales this by the sector's consumption share
	# and baseline aggregate labor income.
	autonomous_demand::Vector{Float64}
	# Investment demand financed by debt (Venue 1 layer). Same form as
	# autonomous_demand; the gap to wage income is the financing record.
	investment_shock::Vector{Float64}
	function Shocks(supply_shock::Vector{Float64}, demand_shock::Vector{Float64},
			demand_shock_raw::Vector{Float64}, autonomous_demand::Vector{Float64},
			investment_shock::Vector{Float64})
		n = length(supply_shock)
		all(length(v) == n for v in
			(demand_shock, demand_shock_raw, autonomous_demand, investment_shock)) ||
			throw(DimensionMismatch("all shock vectors must have the same length"))
		new(supply_shock, demand_shock, demand_shock_raw,
			autonomous_demand, investment_shock)
	end
end

function _shock_vectors(vectors...)
	n = length(first(vectors))
	all(length(v) == n for v in vectors) ||
		throw(DimensionMismatch("all shock vectors must have the same length"))
	map(v -> Float64.(v), vectors)
end

# Backward-compatible positional constructors, including the original 3-arg form.
function Shocks(supply_shock::AbstractVector, demand_shock::AbstractVector,
		demand_shock_raw::AbstractVector, autonomous_demand::AbstractVector,
		investment_shock::AbstractVector)
	a, b, c, d, e = _shock_vectors(supply_shock, demand_shock, demand_shock_raw,
		autonomous_demand, investment_shock)
	Shocks(a, b, c, d, e)
end

Shocks(supply_shock::AbstractVector, demand_shock::AbstractVector,
		demand_shock_raw::AbstractVector) =
	Shocks(supply_shock, demand_shock, demand_shock_raw,
		zeros(length(supply_shock)), zeros(length(supply_shock)))

"""Keyword convenience constructor for optional extra-demand shock vectors."""
function Shocks(supply_shock::AbstractVector, demand_shock::AbstractVector;
		demand_shock_raw=zeros(length(supply_shock)), autonomous_demand=zeros(length(supply_shock)),
		investment_shock=zeros(length(supply_shock)))
	Shocks(supply_shock, demand_shock, demand_shock_raw, autonomous_demand, investment_shock)
end


mutable struct Model{T <: ModelType}
	data::Data
	shocks::Shocks
	options::T
	# Financing closure (Foundation II). Default: NoFinancing (baseline
	# reference). Concrete types live in src/closures/financing/financing.jl.
	financing::AbstractFinancing
end

# 3-arg convenience constructor: defaults to the unfinanced baseline
# reference. Financing experiments construct the 4-arg form.
function Model(data::Data, shocks::Shocks, options)
	Model(data, shocks, options, NoFinancing())
end

labor_closure(model::Model) = labor_closure(model.options)

"""
	assemble_data(io, d) -> Data

Assemble a `Data` object from a raw IO table and the NamedTuple returned by
`generate_data`. The v3 open-economy fields start at their closed-absorption
defaults (zeros / saving rate 0) and are filled by `recalibrate_open` in
`src/core/calibration.jl`; `household_baseline` keeps the legacy Törnqvist
default `consumption_share .* Σ labor_share` so uncalibrated datasets behave
as before (ADR-0005), instead of cbase2's zeros default. Shared by `read_data`
(full table) and `retained_dataset` (rebuilt retained-sector table).
"""
function assemble_data(io::DataFrames.DataFrame, d::NamedTuple)
	n = length(d.λ)
	household_baseline = d.consumption_share .* sum(d.labor_share)
	# A-bill (LaForge fix): the domestic intermediate bill and the intermediate
	# imports, read from the raw table and scaled to model units.
	A_bill, M_int = _domestic_bills(io, n, d.gdp_production)
	return Data(io, d.Ω, d.Ω_raw, d.consumption_share, d.factor_share, d.λ,
			d.labor_share, d.consumption_share_gross_output, d.grossy, d.value_added,
			d.gross_output_basic, d.value_added_components, d.imports_intermediate,
			d.import_share, d.domestic_final_demand, zeros(n), household_baseline,
			zeros(n), zeros(n), zeros(n), 0.0, A_bill, M_int, d.gdp_production,
			d.gdp_income, d.gdp_expenditure)
end

"""
	read_data(filename::String)

Given a filename of a IO table located in the /data directory this returns the CESData, where shocks are set to ones
 given a filename of a io table located in the /data directory this returns the cesdata, where shocks are set to ones
and elasticities are set to the ones presente in the paper by b&f
"""
function read_data(filename::String; datadir::AbstractString = pwd())::Data
	filedir = joinpath(datadir, "data", filename)
	io = CSV.read(filedir, DataFrames.DataFrame, delim = ";", decimal = ',', missingstring = ["-", "x"]) #read in from csv
	DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren) #name the indices after the sectors
	io.Sektoren = replace.(io.Sektoren, r"^\s+" => "") #remove unneccasary whitespaces
	io = coalesce.(io, 0) #set nans to 0

	d = generate_data(io)
	# Assemble the Data object from the §4.1 accounting-consistent transformation
	# (compatibility defaults: see `assemble_data`).
	return assemble_data(io, d)
end
