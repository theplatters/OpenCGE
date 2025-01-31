module BeyondHulten

using NonlinearSolve: NonlinearSolve
using CSV: CSV
using DataFrames
using LineSearches: LineSearches
using LinearAlgebra
using StatsBase
using GLMakie
using XLSX

export CESElasticities, LeontiefElasticies, CobbDouglasElasticities, LeontiefElasticiesLabor, Shocks, Data, Model
export read_data, calculate_investment!
export gross_increase, nominal_increase, real_gdp, nominal_gdp, full_labor_slack
export CBData, read_data_cb, cobb_douglas_costfun, cobb_douglas_consumption, cobb_douglas_intermediary_demand, cobb_douglas_prices!
export solve
export CES, Leontief, CobbDouglas
export gdp

export elasticities_gradient, plot_elasticities, plot_prices, standard_shock, standard_tech_shock, ElasticityGradientSolution, plot_wages
include("leontief.jl")
include("util.jl")

end