struct Solution
	prices::Vector{Float64}
	prices_raw::Vector{Float64}
	quantities::Vector{Float64}
	wages::Vector{Float64}
	wages_raw::Vector{Float64}
	consumption::Vector{Float64}
	numeraire::Float64
	real_gdp::Float64
	nominal_gdp::Float64
	model::Model
end

function Solution(prices_raw, quantities, wages, consumption, numeraire, real_gdp, nominal_gdp, model)
	return Solution(
		prices_raw ./ numeraire,
		prices_raw,
		quantities,
		wages ./ numeraire,
		wages,
		consumption,
		numeraire,
		real_gdp,
		nominal_gdp,
		model)
end

function real_gdp(sol::Solution)::Float64
	return sol.real_gdp
end

function nominal_gdp(sol::Solution)::Float64
	return sol.nominal_gdp
end

function wages(sol::Solution)::Vector{Float64}
	return sol.wages
end

function consumption(sol::Solution)::Vector{Float64}
	sol.consumption
end

function Base.getindex(sol::Solution, ::Colon, sector::String)
	index = findfirst(==(sector), sol.model.data.io.Sektoren)
	return Dict(:prices => sol.prices[index], :quantities => sol.quantities[index], :wages => sol.wages[index], :consumption => sol.consumption[index])
end

struct SectorData
	name::String
    price::Float64
    quantity::Float64
    wage::Float64
    consumption::Float64
end

function eachsector(sol::Solution)
    sectors = sol.model.data.io.Sektoren
    return (SectorData(sectors[i],sol.prices[i], sol.quantities[i], sol.wages[i], sol.consumption[i]) for i in 1:length(sol.prices))
end

function Base.getindex(sol::Solution, sector_index::Int)
    sectors = sol.model.data.io.Sektoren
    return SectorData(
        sectors[sector_index],
        sol.prices[sector_index],
        sol.quantities[sector_index],
        sol.wages[sector_index],
        sol.consumption[sector_index]
    )
end
function Base.getindex(sol::Solution, sector_indices::UnitRange{Int})
    sectors = sol.model.data.io.Sektoren
    return [SectorData(
        sectors[i],
        sol.prices[i],
        sol.quantities[i],
        sol.wages[i],
        sol.consumption[i]
    ) for i in sector_indices]
end
function Base.getindex(sol::Solution, indices::Vector{Int})
    sectors = sol.model.data.io.Sektoren
    return [SectorData(
        sectors[i],
        sol.prices[i],
        sol.quantities[i],
        sol.wages[i],
        sol.consumption[i]
    ) for i in indices]
end