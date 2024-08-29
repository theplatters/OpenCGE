include("cobbdouglas.jl")
"""
	solve(model::Model{LeontiefElasticies}; init)

solves the leontief model
"""
function solve(
	model::Model{LeontiefElasticies};
	init = vcat(ones(length(model.data.grossy)), model.data.λ))
	(; data) = model


	# A = (1 - alpha) Omega^T
	#Eigentlich müssen wir Omega ja jetzt noch erweitern
	q = inv(I - diagm(1 .- data.factor_share) * data.Ω)' * (model.shocks.demand_shock)
	p = ones(length(data.λ))
	df = DataFrames.DataFrame(
		Dict("prices" => p,
			"quantities" => q,
			"sectors" => data.io.Sektoren[1:71],
			"value_added" => model.shocks.demand_shock .* data.consumption_share,
			"value_added_nominal_absolute" => data.grossy .* model.shocks.demand_shock .* data.consumption_share,
			"value_added_nominal_relative" => model.shocks.demand_shock .* data.consumption_share,
			"value_added_absolute" => model.shocks.demand_shock .* data.consumption_share .* data.grossy,
			"value_added_relative" => model.shocks.demand_shock .* data.consumption_share,
		))

	df

end

function solve(
	model::Model{LeontiefElasticiesLabor};
	init = vcat(ones(length(model.data.grossy)), model.data.λ))
	(; data, shocks) = model


	consumption = data.io[1:length(data.consumption_share), 75] ./ sum(data.io[78,2:72])
	wages = (Vector(data.io[78, 2:72]) ./ data.grossy)

	A = vcat( hcat(Matrix(data.io[1:71, 2:72]) ./  (data.grossy'),consumption),
		hcat(wages', 0))

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	shock =  (shocks.demand_shock .- 1) .* consumption
	@info shock
	q = inv(I - A) * (vcat(shock, 0))
	p = ones(length(q))
	df = DataFrames.DataFrame(
		Dict("prices" => p,
			"quantities" => q,
			"sectors" => vcat(data.io.Sektoren[1:71], data.io.Sektoren[78]),
		))

	df

end

"Calculates the gdp of a leontief solution"
gdp(solution, model) = 1 + solution.quantities[72] / sum(model.data.io[findfirst(==("Bruttolöhne und -gehälter"),model.data.io.Sektoren),2:72])
