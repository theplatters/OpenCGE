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
		prices_raw, quantities,
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
