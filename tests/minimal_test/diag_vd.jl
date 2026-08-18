# Diagnostic: run the actual variance_decomposition and pilot on real data.
using CSV: CSV
using DataFrames
using LinearAlgebra
using NonlinearSolve: NonlinearSolve
using LineSearches: LineSearches
using StatsBase
using Printf
using Statistics
using ProgressMeter
using DelimitedFiles

module GLMakie end
module XLSX end

module BeyondHulten
using NonlinearSolve: NonlinearSolve
using CSV: CSV
using DataFrames
using LineSearches: LineSearches
using LinearAlgebra
using StatsBase
using Printf
using Statistics
using ProgressMeter
using DelimitedFiles
export CESElasticities, MobileLaborCESElasticities, Solution, Shocks, Data, Model
export read_data, standard_shock, solve, CES, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, real_gdp, nominal_gdp
export variance_decomposition, summary_table, pilot_eta_sweep
const inflator = 1.46
include("interface.jl")
include("solution.jl")
include("ces.jl")
include("mobile_labor.jl")
include("variance_decomposition.jl")
include("util.jl")
include("impulses.jl")
end

using .BeyondHulten
cd("/workspace/BFrep/(3)BeyondHulten")

data = Data("I-O_DE2019_formatiert.csv")
shocks = standard_shock(data)

println("=== variance_decomposition (default grid) ===")
vd = variance_decomposition(data, shocks;
    η_values = [0.0, 0.5, 1.0, 2.0, 10.0],
    ϵ_values = [0.1, 0.5, 0.99],
    θ_values = [0.1, 0.5, 0.99],
    σ_values = [0.1, 0.5, 0.99],
    output = :real_gdp, verbose=false)
summary_table(vd)

println("\n=== pilot_eta_sweep (reduced grid) ===")
pilot_eta_sweep(data, shocks; θ=0.5, ϵ=0.5, σ=0.9)
