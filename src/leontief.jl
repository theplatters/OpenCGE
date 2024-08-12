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
	q = inv(I - diagm(1 .- data.factor_share) * data.Ω)' * (model.shocks.demand_shock .* data.consumption_share) 
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
