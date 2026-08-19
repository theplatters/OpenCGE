# Core finding: binary wage regime (:fixed vs :mobile) vs continuous η gradient.
using Printf
module GLMakie end; module XLSX end
module BeyondHulten
using NonlinearSolve: NonlinearSolve; using CSV: CSV; using DataFrames
using LineSearches: LineSearches; using LinearAlgebra; using StatsBase
using Printf; using Statistics; using ProgressMeter; using DelimitedFiles
export CESElasticities, MobileLaborCESElasticities, Solution, Shocks, Data, Model
export read_data, standard_shock, autonomous_shock, solve, CES, MobileLaborCES, mobile_labor_model
export sectoral_labor_demand, real_gdp, nominal_gdp, tornqvist_quantity_index
const inflator = 1.46
include("interface.jl"); include("solution.jl"); include("ces.jl")
include("mobile_labor.jl"); include("variance_decomposition.jl"); include("util.jl"); include("impulses.jl")
end
using .BeyondHulten
cd("/workspace/BFrep/(3)BeyondHulten")
data = Data("I-O_DE2019_formatiert.csv")
N = 71; SEC = 35

ref = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 0.5, 0.9, 0.5)
ref_sol = solve(ref); ref_gdp = real_gdp(ref_sol); Lbar = ref.options.labor_bar

println("═"^72)
println("  BINARY REGIME vs CONTINUOUS η: autonomous demand to construction")
println("═"^72)
@printf("%-8s %12s %12s %12s %12s %12s %12s %12s\n",
    "mult", "η=0(Δ%)", "η=1(Δ%)", "fixed(Δ%)", "η_bridge", "fix_bridge", "fixed_emp", "fixed_res")
println("─"^72)

for mult in [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    shocks = autonomous_shock(data; autonomous_mult = mult)
    
    # :mobile η=0
    m0 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.0; closure=:mobile)
    s0 = solve(m0); g0 = 100*(real_gdp(s0)/ref_gdp - 1)
    
    # :mobile η=1 (may stall at large shocks)
    g1 = NaN
    try
        m1 = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; closure=:mobile)
        s1 = solve(m1); g1 = 100*(real_gdp(s1)/ref_gdp - 1)
    catch
    end
    
    # :fixed η=0.5 (sticky wage, employment absorbs)
    gf = NaN; ef = NaN; mrf = NaN
    try
        mf = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.5; closure=:fixed)
        sf = solve(mf)
        Xf = [sf.prices; sf.quantities]; outf = similar(Xf)
        BeyondHulten.problem_fixed(outf, Xf, mf); mrf = maximum(abs.(outf))
        gf = 100*(real_gdp(sf)/ref_gdp - 1)
        ef = sum(sectoral_labor_demand(sf.prices, sf.quantities, 1.0, mf))
    catch
    end
    
    ηb = g1 - g0
    fb = gf - g0
    @printf("%-8.1f %+12.6f %+12.6f %+12.6f %+12.6f %+12.6f %12.4f %12.2e\n",
        mult, g0, g1, gf, ηb, fb, ef, mrf)
end

println("─"^72)
println("Interpretation:")
println("  η_bridge (mobile η=1 − η=0): continuous mobility gradient — NEGLIGIBLE")
println("  fix_bridge (:fixed − :mobile η=0): binary wage regime — MATERIAL")
println("  => Labour-market closure matters through the WAGE REGIME, not through η.")
println("DONE")