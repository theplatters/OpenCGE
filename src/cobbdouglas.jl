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
  α .* p .* y .* data.labor_share .^ -1
end

function cobb_douglas_intermediary_demand(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  r = p .^ data.Ω

  (data.Ω') * (β .* y .* cobb_douglas_costfun(p, y, model)) .* inv.(p)
end

function cobb_douglas_costfun(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  r = p .^ (data.Ω)
  inv.(supply_shock) .* (w .^ α) .* (prod(r, dims=2) .^ β) .* α .^ -α .* prod((β .* data.Ω) .^ (-β .* data.Ω), dims=2)
end

function cobb_douglas_consumption(p, y, model)
  (; data, options, shocks) = model
  (; α, β) = (options.elasticities)
  (; supply_shock, demand_shock) = shocks

  w = cobb_douglas_wages(p, y, model)
  C = w' * data.labor_share
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
  grossy = Vector(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])
  @info grossy
  df = DataFrames.DataFrame(
    Dict("prices" => p,
      "quantities" => q,
      "sectors" => data.io.Sektoren[1:71],
      "value_added_nominal_absolute" => p .* data.labor_share .* cobb_douglas_wages(p, q, model) .* grossy,
      "value_added_nominal_relative" => p .* data.labor_share .* cobb_douglas_wages(p, q, model),
      "value_added_absolute" => data.labor_share .* cobb_douglas_wages(p, q, model) .* grossy,
      "value_added_relative" => data.labor_share .* cobb_douglas_wages(p, q, model),
    ))

  df
end

