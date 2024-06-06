struct CBELasticities <: Elasticities
	α::Float64
	β::Float64
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
	p = @view x[1:N]
	y = @view x[N+1:end]


	out = zeros(eltype(x), 2 * N)

	price!(view(out,1:N), p, y, data)
	out[N+1:end] .= y - intermediary_demand(p, y, data) - consumption(p, y, data)
	out
end

function wages(p, y, data)
	(;α, β) = data.elasticities
	(;supply_shock, demand_shock) = data.shocks
	α * p .* y .* data.labor_share .^ -1
end

function cobb_douglas_intermediary_demand(p, y, data)
	(;α, β) = data.elasticities
	(;supply_shock, demand_shock) = data.shocks

	w = wages(p, y, data)
	r = data.Ω * p


	(β / (α + β)) .* cobb_douglas_costfun(p, y, data) .* r .^ -1 .* (data.Ω * y)
end


function cobb_douglas_costfun(p, y, data)
	(;supply_shock, demand_shock) = data.shocks
	(;α, β) = data.elasticities
	w = wages(p, y, data) 
	r =  (data.Ω * p)
	(α + β) * supply_shock .^ -inv(α + β) .* y .^ inv(α + β)  * (α^(-α/(α + β)) * β^(-β/(α + β))) .* w .^ (α / (α + β)) .* r .^ (β / (α + β))
end

function cobb_douglas_prices!(out, p, y, data::CBData)
	out .= p .- cobb_douglas_costfun(p, y, data)
	nothing
end


function cobb_douglas_consumption(p, y, data)
	(;α, β) = data.elasticities
	(;supply_shock, demand_shock) = data.shocks
	w = wages(p, y, data)
	r = data.Ω * p
	C = w' * data.labor_share
	C * demand_shock .* p .^ (-0.95) .* data.consumption_share
end

function read_data_cb(filename::String)
	filedir = joinpath(pwd(), "data/", filename)
	io = CSV.read(filedir, DataFrames.DataFrame, delim = ";", decimal = ',', missingstring = ["-", "x"]) #Read in from CSV
	DataFrames.rename!(io, Symbol(names(io)[1]) => :Sektoren) #Name the indices after the sectors
	io.Sektoren = replace.(io.Sektoren, r"^\s+" => "") #Remove unneccasary whitespaces
	io = coalesce.(io, 0) #Set NANS to 0

	Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, grossy = generateData(io)
	#return a mutable structure element (see above):
	return CBData(io, Ω, consumption_share, factor_share, λ, labor_share, consumption_share_go, Shocks(ones(71), ones(71)), CBELasticities(0.3,0.7), grossy)
end

function solve_cobb_douglas_modell(
	data::Data,
	shocks::Shocks,
    elasticities::Elasticities,
	init = Complex.(vcat(ones(length(data.λ)), data.λ)))

    set_shocks!(data, shocks)
    set_elasticities!(data, elasticities)

	f = NonlinearSolve.NonlinearFunction((x, u) -> generalized_problem(x, u, cobb_douglas_prices!, cobb_douglas_intermediary_demand, cobb_douglas_consumption))
	prob = NonlinearSolve.NonlinearProblem(f, init, data)

	sol = NonlinearSolve.solve(prob)
end
