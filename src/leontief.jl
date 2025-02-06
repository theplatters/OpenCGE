"""
	solve(model::Model{LeontiefElasticies}; init)

solves the leontief model
"""
function solve(model::Model{Leontief})

  (; data, shocks) = model
  consumption_share = data.io[1:length(data.consumption_share), 75] ./ sum(data.io[78, 2:72])
  consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
                sum |>
                x -> getindex(x, 1:71)

  shock = (shocks.demand_shock .- 1) .* consumption

  if model.options.labor_effect
    wages = (Vector(data.io[78, 2:72]) ./ data.grossy)

    A = vcat(hcat(Matrix(data.io[1:71, 2:72]) ./ (data.grossy'), consumption_share),
      hcat(wages', 0))


    @info shock
    q = inv(I - A) * (vcat(shock, 0))
    p = ones(length(q))


  else

    # A = (1 - alpha) Omega^T
    #Eigentlich müssen wir Omega ja jetzt noch erweitern
    q = inv(I - diagm(1 .- data.factor_share) * data.Ω)' * (shock)
    p = ones(length(q))

  end

  df = DataFrames.DataFrame(
    Dict("prices" => p,
      "quantities" => q,
      "sectors" => model.options.labor_effect ? vcat(data.io.Sektoren[1:71], data.io.Sektoren[78]) : data.io.Sektoren[1:71],
    ))

  return df
end

#=
function solve(
	model::Model{Leontief};
	init = vcat(ones(length(model.data.grossy)), model.data.λ))


	consumption = data.io[1:length(data.consumption_share), 75] ./ sum(data.io[78, 2:72])
	wages = (Vector(data.io[78, 2:72]) ./ data.grossy)

	A = vcat(hcat(Matrix(data.io[1:71, 2:72]) ./ (data.grossy'), consumption),
		hcat(wages', 0))

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	shock = (shocks.demand_shock .- 1) .* consumption
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
=#
"Calculates the gdp of a leontief solution"
gdp(solution, model::Model{Leontief}) = 1 + solution.quantities[72] / sum(model.data.io[findfirst(==("Bruttolöhne und -gehälter"), model.data.io.Sektoren), 2:72])
