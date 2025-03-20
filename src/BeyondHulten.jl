module BeyondHulten

using NonlinearSolve: NonlinearSolve
using CSV: CSV
using DataFrames
using LineSearches: LineSearches
using LinearAlgebra
using StatsBase
using GLMakie
using XLSX

# Export Elasticity Types
export CESElasticities, LeontiefElasticies, CobbDouglasElasticities, LeontiefElasticiesLabor

# Export Data, Models, and Shocks
export Solution
export Shocks, Data, Model, CBData, read_data_cb

# Export Data Reading and Investment Functions
export read_data, calculate_investment!
export gross_increase, nominal_increase, real_gdp, nominal_gdp, full_labor_slack, gdp

# Export Cobb–Douglas Specific Functions
export cobb_douglas_costfun, cobb_douglas_consumption, cobb_douglas_intermediary_demand, cobb_douglas_prices!

# Export Model Constructors & Solver
export solve
export CES, Leontief, CobbDouglas, Solution


# Export Visualization and Utility Functions
export elasticities_gradient, plot_elasticities, standard_shock, standard_tech_shock, ElasticityGradientSolution, plot_wages, impulse_shock, plot_consumption, load_impulses

include("interface.jl")
include("solution.jl")
include("cobbdouglas.jl")
include("leontief.jl")
include("ces.jl")
include("util.jl")
include("impulses.jl")

end
