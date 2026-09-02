# Milestone E (bridge) diagnostic: repurpose eta as B&F (2019) labor-reallocation
# parameter beta. Compare immobile (eta=0) vs mobile (eta=1) under a construction
# (sector 35) productivity (supply) shock. Confirm the mobility bridge appears as a
# real-GDP difference, and that the package residuals stay near zero.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
N = 71
const SEC = 35   # construction sector

function maxres(model, sol)
    maximum(abs, equilibrium_residuals(sol))
end

# B&F bridge case: sector-specific POSITIVE PRODUCTIVITY (supply) shock.
# Use epsilon = 2.0 (>1, substitution regime) so the reallocation gain has the
# correct (positive) sign: immobile labor is costly, mobile labor reallocates.
ss = ones(N); ss[SEC] *= 1.30
shocks = Shocks(ss, ones(N), zeros(N))
base_model = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 2.0, 0.9, 0.0)
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)
@printf("baseline real_gdp = %.6f\n", base_gdp)

etas = [0.0, 0.25, 0.5, 0.75, 1.0]
gdps = Float64[]
resids = Float64[]
wages = Float64[]
sumLs = Float64[]
println("eta   maxRes       realGDP     realGDPΔ%    wage      sumL")
for η in etas
    model = mobile_labor_model(data, shocks, 0.5, 2.0, 0.9, η)
    sol = solve(model)
    mr = maxres(model, sol)
    g = real_gdp(sol)
    sl = sum(sectoral_labor_demand(
        sol.prices_raw, sol.quantities, sol.wages_raw[1], model))
    push!(gdps, g); push!(resids, mr); push!(wages, sol.wages[1]); push!(sumLs, sl)
    @printf("%5.2f %.2e  %10.6f  %+.6f%%  %.6f  %.6f\n",
        η, mr, g, 100*(g/base_gdp - 1), sol.wages[1], sl)
end
g_imm = gdps[1]; g_mob = gdps[end]
@printf("BRIDGE (mobile - immobile real-GDP Δ%%) = %+.6f%%  (≈ %.4f pp)\n",
    100*(g_mob/base_gdp - g_imm/base_gdp), 100*(g_mob/base_gdp - g_imm/base_gdp))
println("DONE")
