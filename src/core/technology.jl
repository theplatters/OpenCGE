# src/core/technology.jl — CES + Leontief + Cobb-Douglas technology (concatenated verbatim:
# src/ces.jl + src/leontief.jl + src/cobbdouglas.jl, in that order)

"""
	calculate_investment!(shock::Shocks, data::Data, investment::Number, sector::String)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Vector{<:Number}, sector)

	consumption = eachcol(data.io[:, DataFrames.Between("Konsumausgaben der privaten Haushalte im Inland", "Exporte")]) |>
				  sum |>
				  x -> getindex(x, 1:71)

	for i in eachindex(sector)
		sector_number = findfirst(==(sector[i]), data.io.Sektoren)
		shocks.demand_shock[sector_number] = 1 + investment[i] / consumption[sector_number]
		println("Demand shock to sector $(sector[i]): $(shocks.demand_shock[sector_number])")
	end
end

"""
	calculate_investment!(shock::Shocks, data::Data, investment::Dict)

Alters the shock vector, so that the shock in the given sector reflects the investment in thousend €
"""
function calculate_investment!(shocks::Shocks, data::AbstractData, investment::Dict{String, Number})
	for (sector, investment) in investment
		calculate_investment!(shocks, data, [investment], [sector])
	end
end

"""
	full_labor_slack_alt(model::Model)

Returns the labor vector adjusted, so that labor can be freely reallocated to accomodate for demand shocks
"""
function full_labor_slack_alt(model::Model)
	(; data, shocks) = model

	q = inv(I - diagm(1 .- data.factor_share) * data.Ω_raw)' * shocks.demand_shock_raw
	data.labor_share + (q ./ sum(data.value_added)) .* (Vector(data.io[findfirst(==("Arbeitnehmerentgelt im Inland"), data.io.Sektoren), 2:72]) ./ data.grossy)
end

"""
	full_labor_slack(model::Model)

Returns the labor vector adjusted, so that labor can be freely reallocated to accomodate for demand shocks
"""
function full_labor_slack(model::Model)
	(; data, shocks) = model

	q = inv(I - diagm(1 .- data.factor_share) * data.Ω_raw)' * shocks.demand_shock_raw
	data.labor_share + inv(I - diagm(1 .- data.factor_share) * data.Ω_raw)' * (data.consumption_share_gross_output .* ((shocks.demand_shock .* data.labor_share) - data.labor_share))
end


function empirical_labor_slack(model::Model, unemployment_rate::Float64 = 0.031)
	(; data, shocks) = model
	(1 / (1 - unemployment_rate)) * data.labor_share
end

"""
  problem(X, model::Model{CES})

The objective function as specified in B&F with the added demand shocks, X is the 2*sectors-sized vector,
data contains the parameters and labor_reallocation is a function that specifies how labor is reallocated accross sectors
"""
function problem(out::Vector, X::Vector, model::Model{CES})

	(; data, options, shocks) = model
	N = length(data.factor_share)
	p = max.(X[1:N], 0)
	y = max.(X[N+1:end], 0)

	(; supply_shock, demand_shock) = shocks
	(; consumption_share, Ω_raw, factor_share) = data
	(; ϵ, θ, σ) = options.elasticities
	labor = options.labor_slack(model)



	intermediate_price = (Ω_raw * p .^ (1 - θ)) .^ (1 / (1 - θ))

	cpi = sum(data.consumption_share .* p .^ (1 - σ))^(1 / (1 - σ))
	w = p .* (supply_shock .^ ((ϵ - 1) / ϵ)) .* (factor_share .^ (1 / ϵ)) .* (y .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)

	C = w' * labor
	final_demand = (C * p .^ (-σ) .* demand_shock .* consumption_share) ./ cpi .^ (-σ)
	intermediary_demand = p .^ (-θ) .* (Ω_raw' * (p .^ ϵ .* supply_shock .^ (ϵ - 1) .* intermediate_price .^ (θ - ϵ) .* (1 .- factor_share) .* y))
	out[1:N] .= p - (supply_shock .^ (ϵ - 1) .* (factor_share .* w .^ (1 - ϵ) + (1 .- factor_share) .* intermediate_price .^ (1 - ϵ))) .^ (1 / (1 - ϵ))
	out[N+1:end] .= y - intermediary_demand - final_demand
	nothing
end

"""Return the exact residual vector for a legacy CES equilibrium."""
function equilibrium_residuals(model::Model{CES}, X::AbstractVector)
	N = length(model.data.factor_share)
	length(X) == 2N || throw(DimensionMismatch("CES expects a 2N vector"))
	out = zeros(Float64, 2N)
	problem(out, collect(X), model)
	out
end

# `Solution` is defined in src/core/equilibrium.jl (included after this file),
# so the method below duck-types its second argument; behavior is unchanged.
_equilibrium_residuals(::Model{CES}, sol) =
	equilibrium_residuals(sol.model, [sol.prices_raw; sol.quantities])

"""
	solve_ces_model(model::Model{CES}; init = [ones(length(model.data.λ)); model.data.λ])
The main function of this module, input the relavant model data, shocks and optionally labor_reallocation and
starting vectors and get back the simulated adapted prices and quantities

"""
function solve(
	model::Model{CES};
	init = [ones(length(model.data.λ)); model.data.λ],
)
	(; data, options, shocks) = model


	#defines the function:
	#defines the concrete problem to be solved (i.e. with inserted parameter values):
	ProbN = NonlinearSolve.NonlinearProblem(problem, init, model)
	x = NonlinearSolve.solve(ProbN, reltol = 1e-8, abstol = 1e-8).u


	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]


	labor = options.labor_slack(model)
	(; ϵ, θ, σ) = options.elasticities
	wages = p .* (shocks.supply_shock .^ ((ϵ - 1) / ϵ)) .* (data.factor_share .^ (1 / ϵ)) .* (q .^ (1 / ϵ)) .* labor .^ (-1 / ϵ)
	consumption_share = shocks.demand_shock .* data.consumption_share

	numeraire = (data.consumption_share' * p .^ (1 - σ))^(1 / (1 - σ))
	consumption = wages' * labor .* consumption_share .* (p / numeraire) .^ (-σ)
	real_gdp = tornqvist_quantity_index(
		p,
		consumption,
		ones(length(p)),
		data.consumption_share,
	)

	nominal_gdp = wages' * labor
	return Solution(p, q, wages, consumption, numeraire, real_gdp, nominal_gdp, model)
end

# --- src/leontief.jl (verbatim) ---

"""
	solve(model::Model{LeontiefElasticies}; init)

solves the leontief model
"""
function solve(model::Model{Leontief})

	(; data, shocks) = model
	consumption_share = data.io[1:length(data.consumption_share), 75] ./ sum(data.io[78, 2:73])

	shock =  shocks.demand_shock_raw

	wages = (Vector(data.io[78, 2:72]) ./ data.grossy)

	A = vcat(hcat(Matrix(data.io[1:71, 2:72]) ./ (data.grossy'), consumption_share),
		hcat(wages', 0))


	q = inv(I - A) * (vcat(shock, 0))
	p = ones(length(q))


	value_added = Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])
	value_added_share = value_added ./ Vector(data.io[findfirst(==("Produktionswert"), data.io.Sektoren), 2:72])
	real_gdp =
		1 +
		sum(value_added_share .* q[1:71]) ./
		sum(value_added)
	q = [data.λ;0] .+  q ./ sum(value_added)
	return Solution(p, q, ones(length(q)), shocks.demand_shock + q[72] .* consumption_share, 1, real_gdp, real_gdp, model)
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

# --- src/cobbdouglas.jl (verbatim) ---

function generalized_problem(x, model, costfun, intermediary_demand, consumption)

  N = length(model.data.λ)
  p = max.(0, x[1:N])
  y = max.(0, x[N+1:end])


  out = zeros(eltype(x), 2 * N)

  out[1:N] .= p .- costfun(p, y, model)
  out[N+1:end] .= y - intermediary_demand(p, y, model) - consumption(p, y, model)
  out
end

function cobb_douglas_wages(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  labor  = options.labor_slack(model)
  α .* p .* y .* labor .^ -1
end

function cobb_douglas_intermediary_demand(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  r = p .^ data.Ω_raw

  (data.Ω_raw') * (β .* y .* cobb_douglas_costfun(p, y, model)) .* inv.(p)
end

function cobb_douglas_costfun(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  r = p .^ (data.Ω_raw)
  inv.(supply_shock) .* (w .^ α) .* (prod(r, dims=2) .^ β) .* α .^ -α .* prod((β .* data.Ω_raw) .^ (-β .* data.Ω_raw), dims=2)
end

function cobb_douglas_consumption(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  labor  = options.labor_slack(model)
  C = w' * labor 
  C * demand_shock .* p .^ (-1) .* data.consumption_share
end

function solve(
  model::Model{CobbDouglas};
  init=(vcat(ones(length(model.data.λ)), model.data.λ)))

  (; data) = model

  f = NonlinearSolve.NonlinearFunction((x, u) -> generalized_problem(x, u, cobb_douglas_costfun, cobb_douglas_intermediary_demand, cobb_douglas_consumption))
  prob = NonlinearSolve.NonlinearProblem(f, init, model)

  x = NonlinearSolve.solve(prob)
  p = x[1:length(data.consumption_share)]
  q = x[(length(data.consumption_share)+1):end]
  wages = cobb_douglas_wages(p,q,model)
  labor = model.options.labor_slack(model)
  consumption = cobb_douglas_consumption(p,q,model)
  real_gdp = tornqvist_quantity_index(
    p,
    consumption,
    ones(length(p)),
    data.consumption_share,
  )
  numeraire = mean(p, weights(consumption))
  grossy = Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])

  return Solution(p, q, wages, consumption, numeraire, real_gdp, wages' * labor, model)

end
