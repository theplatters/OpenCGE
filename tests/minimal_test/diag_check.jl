using Printf
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
export read_data, standard_shock, autonomous_shock, solve, CES, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, real_gdp, nominal_gdp, tornqvist_quantity_index
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
sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"
ci = findfirst(==(sector), data.io.Sektoren)
println("construction index = ", ci)
println("nsectors = ", length(data.consumption_share))
sh = autonomous_shock(data; autonomous_mult = 0.2)
println("autonomous_demand nonzero indices: ", findall(!iszero, sh.autonomous_demand))
println("autonomous_demand[ci] = ", sh.autonomous_demand[ci])
println("consumption_share[ci] = ", data.consumption_share[ci])
println("sum(labor_share) = ", sum(data.labor_share))
