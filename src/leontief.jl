include("cobbdouglas.jl")

function leontief_intermediary_demand(p, y, data)
	return (diagm(1 .- data.factor_share) * data.Ω)' * y
end

function leontief_costfun(p, y, data)
	return data.Ω * p
end

function leontief_consumption(p, y, data)
	return data.shocks.demand_shock .* data.consumption_share
end

function solve_leontief_modell(data,
	shocks::Shocks,
	init = vcat(ones(length(data.grossy)), data.λ))

	set_shocks!(data, shocks)

	f = NonlinearSolve.NonlinearFunction((x, u) -> generalized_problem(x, u, leontief_costfun, leontief_intermediary_demand, leontief_consumption))
	prob = NonlinearSolve.NonlinearProblem(f, init, data)

	x = NonlinearSolve.solve(prob)
	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]
	df = DataFrames.DataFrame(
		Dict("prices" => p,
			"quantities" => q,
			"sectors" => data.io.Sektoren[1:71],
			"value_added" => data.labor_share .* cobb_douglas_wages(p, q, data)),
	)

	df

end