abstract type AbstractElasticities end
abstract type AbstractData end


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
CES(elasticities::CESElasticities) = CES(elasticities, full_labor_slack, false)
CES(elasticities::CESElasticities, labor_slack) = CES(elasticities, labor_slack, false)
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
	# Conditional intermediate-input shares in B&F orientation:
	# rows are using sectors and columns are supplying sectors.
	Ω::Matrix{Float64}
	consumption_share::Vector{Float64}
	factor_share::Vector{Float64}
	λ::Vector{Float64}
	labor_share::Vector{Float64}
	consumption_share_gross_output::Vector{Float64}
	grossy::Vector{Float64}
	value_added::Vector{Float64}
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

"""
	generateData(io::DataFrames.DataFrame)

Helper function that pulls out the key econometric variables used in the model out of the extended io table
# Example
```julia-repl
julia> Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go = generateData(io)
```
"""
function generate_data(io::DataFrames.DataFrame)
	number_sectors = 71
	# Destatis reports supplied products in rows and using production sectors in
	# columns. B&F's IO matrix uses the opposite orientation before normalizing
	# every using sector's intermediate bundle to one.
	intermediate_use = Matrix(coalesce.(io[1:number_sectors, 2:number_sectors+1], 0.0))
	Ω = conditional_input_shares(intermediate_use)

	grossy = io[1:number_sectors, "Gesamte Verwendung von Gütern"]
	consumption = eachcol(io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:number_sectors)
	value_added = Vector(io[findfirst(==("Bruttowertschöpfung"), io.Sektoren), 2:number_sectors+1])

	factor_share = value_added ./ grossy
	consumption_share = (I - diagm(1 .- factor_share) * Ω)' * grossy
	@views consumption_share[consumption_share.<0] .= 0
	consumption_share = consumption_share / sum(consumption_share)
	λ = (inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share)
	labor_share = λ .* factor_share
	consumption_share_gross_output = consumption ./ grossy
	return Ω, consumption_share, factor_share, λ, labor_share, consumption_share_gross_output, grossy, value_added
end


""" 
Shocks

Two vectors, that have to be the same length as the amount of sectors.
Each entry that differs from one, represents a percentage shock in that sector on demand/supply.

# Example
```julia-repl
julia> Shocks(ones(76),ones(76))
```
"""
struct Shocks
	supply_shock::Vector{Float64}
	demand_shock::Vector{Float64}
	demand_shock_raw::Vector{Float64}
end


mutable struct Model{T <: ModelType}
	data::Data
	shocks::Shocks
	options::T
end

"""
	read_data(filenem::String)

Given a filename of a IO table located in the /data directory this returns the CESData, where shocks are set to ones
 given a filename of a io table located in the /data directory this returns the cesdata, where shocks are set to ones
and elasticities are set to the ones presente in the paper by b&f
"""
function read_data(filename::String)::Data
	filedir = joinpath(pwd(), "data/", filename)
	io = CSV.read(filedir, DataFrames.DataFrame, delim = ";", decimal = ',', missingstring = ["-", "x"]) #read in from csv
	DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren) #name the indices after the sectors
	io.Sektoren = replace.(io.Sektoren, r"^\s+" => "") #remove unneccasary whitespaces
	io = coalesce.(io, 0) #set nans to 0

	Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, grossy, value_added = generate_data(io)
	#return a mutable structure element (see above):	
	return Data(io, Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, grossy, value_added)
end
