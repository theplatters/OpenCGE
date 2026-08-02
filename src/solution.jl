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

"""
	tornqvist_quantity_index(prices, quantities, base_prices, base_quantities)

Compute the discrete Törnqvist approximation to the Divisia quantity index.
The logarithmic growth rate is the expenditure-share-weighted change in final
quantities,

`log(Q₁/Q₀) = sumᵢ 0.5 * (sᵢ₀ + sᵢ₁) * log(qᵢ₁/qᵢ₀)`.

Goods with zero expenditure shares in both equilibria are omitted. A good with
a positive expenditure share in either equilibrium must have a strictly
positive quantity in both equilibria.
"""
function tornqvist_quantity_index(
	prices::AbstractVector{<:Real},
	quantities::AbstractVector{<:Real},
	base_prices::AbstractVector{<:Real},
	base_quantities::AbstractVector{<:Real},
)::Float64
	lengths = length.((prices, quantities, base_prices, base_quantities))
	all(==(first(lengths)), lengths) ||
		throw(DimensionMismatch("prices and quantities must have matching lengths"))
	any(x -> x < 0, prices) && throw(DomainError(prices, "prices must be nonnegative"))
	any(x -> x < 0, base_prices) && throw(DomainError(base_prices, "base prices must be nonnegative"))
	any(x -> x < 0, quantities) && throw(DomainError(quantities, "quantities must be nonnegative"))
	any(x -> x < 0, base_quantities) && throw(DomainError(base_quantities, "base quantities must be nonnegative"))

	current_expenditure = prices .* quantities
	base_expenditure = base_prices .* base_quantities
	current_total = sum(current_expenditure)
	base_total = sum(base_expenditure)
	current_total > 0 || throw(DomainError(current_total, "current expenditure must be positive"))
	base_total > 0 || throw(DomainError(base_total, "base expenditure must be positive"))

	current_shares = current_expenditure ./ current_total
	base_shares = base_expenditure ./ base_total
	active = (current_shares .+ base_shares) .> 0
	any(quantities[active] .<= 0) &&
		throw(DomainError(quantities, "positive-share goods require positive quantities"))
	any(base_quantities[active] .<= 0) &&
		throw(DomainError(base_quantities, "positive-share goods require positive base quantities"))

	log_growth = sum(
		0.5 .* (base_shares[active] .+ current_shares[active]) .*
		log.(quantities[active] ./ base_quantities[active]),
	)
	return exp(log_growth)
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
	return (SectorData(sectors[i], sol.prices[i], sol.quantities[i], sol.wages[i], sol.consumption[i]) for i in 1:length(sol.prices))
end

function Base.getindex(sol::Solution, sector_index::Int)
	sectors = sol.model.data.io.Sektoren
	return SectorData(
		sectors[sector_index],
		sol.prices[sector_index],
		sol.quantities[sector_index],
		sol.wages[sector_index],
		sol.consumption[sector_index],
	)
end
function Base.getindex(sol::Solution, sector_indices::UnitRange{Int})
	sectors = sol.model.data.io.Sektoren
	return [SectorData(
		sectors[i],
		sol.prices[i],
		sol.quantities[i],
		sol.wages[i],
		sol.consumption[i],
	) for i in sector_indices]
end
function Base.getindex(sol::Solution, indices::Vector{Int})
	sectors = sol.model.data.io.Sektoren
	return [SectorData(
		sectors[i],
		sol.prices[i],
		sol.quantities[i],
		sol.wages[i],
		sol.consumption[i],
	) for i in indices]
end

function multiplier(sol::Solution)::Float64
	(; data, shocks) = sol.model
	simple_effect = 1 + sum(shocks.demand_shock_raw) ./ sum(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])
	(1 .- sol.real_gdp) ./ (1 .- simple_effect)
end
