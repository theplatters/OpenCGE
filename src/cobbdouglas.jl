struct CBELasticities <: Elasticities
	α::Vector{Float64}
	β::Vector{Float64}
end


mutable struct CBData <: Data
	io::DataFrames.DataFrame
	Ω::Matrix{Float64}
	consumption_share::Vector{Float64}
	factor_share::Vector{Float64}
	λ::Vector{Float64}
	labor_share::Vector{Float64}
	consumption_share_gross_output::Vector{Float64}
	shocks::Shocks
	elasticities::CBELasticities
	grossy::Vector{Float64}
end

function generalized_problem(x, data::Data, price!, intermediary_demand, consumption)
	N = length(data.λ)
	p = max.(0,x[1:N])
	y = max.(0,x[N+1:end])


	out = zeros(eltype(x), 2 * N)

	price!(view(out, 1:N), p, y, data)
	out[N+1:end] .= y - intermediary_demand(p, y, data) - consumption(p, y, data)
	out
end

function wages(p, y, data)
	(; α, β) = (data.elasticities)
	(; supply_shock, demand_shock) = data.shocks
	α .* p .* y .* data.labor_share .^ -1
end

function cobb_douglas_intermediary_demand(p, y, data)
	(; α, β) = (data.elasticities)
	(; supply_shock, demand_shock) = data.shocks

	w = wages(p, y, data)
	r = p .^ data.Ω

	(data.Ω') * (β .* y .* cobb_douglas_costfun(p, y, data)) .* inv.(p)
end


function cobb_douglas_costfun(p, y, data)
	(; α, β) = (data.elasticities)
	(; supply_shock, demand_shock) = data.shocks
	w = wages(p, y, data)
	r = p .^ data.Ω
	inv.(supply_shock) .* (w .^ α) .* (prod(r, dims = 2) .^ β)
end

function cobb_douglas_prices!(out, p, y, data::CBData)
	out .= (p .- cobb_douglas_costfun(p, y, data))
	nothing
end


function cobb_douglas_consumption(p, y, data)
	(; α, β) = data.elasticities
	(; supply_shock, demand_shock) = data.shocks
	w = wages(p, y, data)
	r = p .^ data.Ω
	C = w' * data.labor_share
	C * demand_shock .* p .^ (-0.95) .* data.consumption_share #taking a CES consumption function for now
end

function read_data_cb(filename::String)
	filedir = joinpath(pwd(), "data/", filename)
	io = CSV.read(filedir, DataFrames.DataFrame, delim = ";", decimal = ',', missingstring = ["-", "x"]) #Read in from CSV
	DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren) #Name the indices after the sectors
	io.Sektoren = replace.(io.Sektoren, r"^\s+" => "") #Remove unneccasary whitespaces
	io = coalesce.(io, 0) #Set NANS to 0

	Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, grossy = generateData(io)
	#return a mutable structure element (see above):
	return CBData(io, Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, Shocks(ones(71), ones(71)), CBELasticities(factor_share, 1 .- factor_share), grossy)
end

function solve_cobb_douglas_modell(
	data::Data,
	shocks::Shocks,
	elasticities::Elasticities,
	init = (vcat(ones(length(data.λ)), data.λ)))

	set_shocks!(data, shocks)
	set_elasticities!(data, elasticities)

	f = NonlinearSolve.NonlinearFunction((x, u) -> generalized_problem(x, u, cobb_douglas_prices!, cobb_douglas_intermediary_demand, cobb_douglas_consumption))
	prob = NonlinearSolve.NonlinearProblem(f, init, data)

	x = NonlinearSolve.solve(prob)
	p = x[1:length(data.consumption_share)]
	q = x[(length(data.consumption_share)+1):end]
	df = DataFrames.DataFrame(
		Dict("prices" => p,
			"quantities" => q,
			"sectors" => data.io.Sektoren[1:71],
			"value_added" => data.labor_share .* wages(p, q, data)),
	)

	df
end
