module BeyondHulten

import NonlinearSolve
import CSV
import DataFrames
import LineSearches
import SciMLNLSolve
using LinearAlgebra

export CESElasticities, LeontiefElasticies, CobbDouglasElasticities, LeontiefElasticiesLabor, Shocks, Data, Model
export read_data, calculate_investment! 
export gross_increase, nominal_increase,real_gdp, nominal_gdp, full_demand_labor_allocation
export CBData, read_data_cb, cobb_douglas_costfun, cobb_douglas_consumption, cobb_douglas_intermediary_demand, cobb_douglas_prices!
export solve
include("leontief.jl")


end