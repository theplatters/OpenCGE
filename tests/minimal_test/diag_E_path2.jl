# Path 2 verification: sticky-wage / unemployment closure.
# Compare :mobile (wage market-clears, wage pinned -> employment fixed) vs
# :fixed (wage sticky at 1 -> employment absorbs the shock).
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")

function maxres(model, sol)
    maximum(abs, equilibrium_residuals(sol))
end

# Baseline (no shock)
base = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.5), data))
base_sol = solve(base)
base_gdp = real_gdp(base_sol)
Lbar = base.options.labor_bar

function run(closure, mult)
    shocks = autonomous_shock(data; autonomous_mult = mult)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 0.5; closure = closure)
    sol = solve(model)
    mr = maxres(model, sol)
    L = sum(sectoral_labor_demand(
        sol.prices_raw, sol.quantities, sol.wages_raw[1], model))
    @printf("%-7s  maxRes=%.2e  realGDPΔ=%+.6f%%  emp=%.6f  unemp=%.6f  wage=%.6f\n",
        closure, mr, 100*(real_gdp(sol)/base_gdp - 1), L, Lbar - L, sol.wages[1])
end

println("Construction autonomous shock (mult=0.2). Baseline emp=", Lbar, " wage=1.0")
print("  :mobile  "); run(:mobile, 0.2)
print("  :fixed   "); run(:fixed, 0.2)
println("\nConstruction autonomous shock (mult=0.5)")
print("  :mobile  "); run(:mobile, 0.5)
print("  :fixed   "); run(:fixed, 0.5)
println("DONE")
