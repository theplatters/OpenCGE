struct CESTemporalElasticities <: AbstractElasticities
	θ::Float64
	ϵ::Float64
	σ::Float64
	intertemporal_elasticity::Float64
end

struct CESTemporalParams
	htm_share::Float64
end

struct CESTemporal <: ModelType
	elasticites::CESTemporalElasticites
	labor_slack::Union{Function, Symbol}
end

struct CESTemporalData <: AbstractData
    Ω::Matrix{Float64}
    consumption_share::Vector{Float64}
    factor_share::Vector{Float64}
    labor_share::Vector{Float64}
end


function problem(out::Vector, X::Vector, model::Model{CESTemporal})
	(; data, options, shocks) = model
	N = length(data.factor_share)
	p = max.(X[1:N], 0)
	y = max.(X[N+1:end], 0)

	(; supply_shock, demand_shock) = shocks
	(; consumption_share, Ω, factor_share) = data
	(; ϵ, θ, σ, intertemporal_elasticity) = options.elasticities
	consumption_share = (demand_shock .^ θ * Ω[1, k] * (p[k] / p[member(1, Ind)])^(1 - theta[member(1, Ind)]) * A[member(1, Ind)]^(theta[member(1, Ind)] - 1))

end
