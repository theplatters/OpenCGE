module BeyondHulten

using NonlinearSolve: NonlinearSolve
using CSV: CSV
using DataFrames
using LineSearches: LineSearches
using LinearAlgebra
using StatsBase
using GLMakie
using XLSX
using DelimitedFiles
# Export Elasticity Types
export CESElasticities, LeontiefElasticies, CobbDouglasElasticities, LeontiefElasticiesLabor, MobileLaborCESElasticities
export AbstractLaborClosure, ExogenousLaborClosure, FlexibleWageClosure, FixedWageClosure, labor_closure

# Export Data, Models, and Shocks
export Solution, eachsector
export Shocks, Data, Model, CBData, read_data_cb

# Export Data Reading and Investment Functions
export read_data, calculate_investment!
export gross_increase, nominal_increase, real_gdp, nominal_gdp, tornqvist_quantity_index, full_labor_slack, gdp

# Export Cobb–Douglas Specific Functions
export cobb_douglas_costfun, cobb_douglas_consumption, cobb_douglas_intermediary_demand, cobb_douglas_prices!

# Export Model Constructors & Solver
export solve
export equilibrium_residuals
export CES, Leontief, CobbDouglas, Solution, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, economy_wide_wage


# Export Visualization and Utility Functions
export elasticities_gradient, standard_shock, standard_tech_shock, autonomous_shock, ElasticityGradientSolution, plot_wages, impulse_shock, plot_consumption, load_impulses, multiplier, full_labor_slack_alt, inflator, cpi
export plot_real_gdp_gradient,plot_nominal_gdp_gradient, panel, diff_lambda, comparison_between_labor_slacks,effect_of_different_elasticities
export eta_sweep, eta_sweep_full, EtaSweepResult, real_gdp_sweep, nominal_gdp_sweep, sectoral_quantities, sectoral_prices
export variance_decomposition, SobolResult, VarianceDecompositionResult, summary_table, eta_sweep_diagnostics, pilot_eta_sweep
const inflator = 1.46

include("interface.jl")
include("solution.jl")
include("cobbdouglas.jl")
include("leontief.jl")
include("ces.jl")
include("mobile_labor.jl")
include("util.jl")
include("impulses.jl")
include("plots.jl")
include("variance_decomposition.jl")

end
