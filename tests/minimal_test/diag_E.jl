# Milestone E verification (Venue 2, then Venue 2+1 layer): an EXPANSIVE
# autonomous export-demand shock should engage the labor channel so that the
# real-GDP response varies with the intersectoral-reallocation parameter η.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")

# Reconstruct the FULL corrected residual vector (incl. autonomous/investment).
function full_res(model, sol)
    N = length(model.data.factor_share)
    p = sol.prices_raw
    residuals = equilibrium_residuals(sol)
    L_i = sectoral_labor_demand(p, sol.quantities, sol.wages_raw[1], model)
    total_income = sol.wages_raw[1] * sum(L_i)
    ax = model.shocks.autonomous_demand
    gv = model.shocks.investment_shock
    cons_base = sum(data.labor_share)
    A = ax .* data.consumption_share .* cons_base; G = gv .* data.consumption_share .* cons_base
    debt = sum(p .* (sol.consumption .+ A .+ G)) - total_income
    return (zp=residuals[1:N], mc=residuals[N+1:2N-1],
        labor=residuals[2N], num=residuals[2N+1], debt=debt)
end

# True no-shock baseline (all shocks = 1 / 0).
base_model = Model(data, Shocks(ones(71), ones(71), zeros(71)), MobileLaborCES(MobileLaborCESElasticities(0.5,0.5,0.9,0.0), data))
base_sol = solve(base_model)
base_gdp = real_gdp(base_sol)

println("=== Milestone E (Venue 2: autonomous export demand) ===")
println("eta     maxRes        realGDPΔ%     wage        debt(financing)")
for η in [0.0, 0.5, 1.0, 5.0, 50.0]
    shocks = autonomous_shock(data; autonomous_mult = 0.2)   # demand_shock=1, autonomous[constr]=0.2
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_res(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f  %.6f\n",
        η, max(maximum(abs.(r.zp)), maximum(abs.(r.mc)), abs(r.labor), abs(r.num)),
        100*(g/base_gdp - 1), sol.wages[1], r.debt)
end

println("\n=== Milestone E (Venue 2 + Venue 1 layer: debt-financed investment) ===")
for η in [0.0, 1.0, 50.0]
    shocks = autonomous_shock(data; investment_mult = 0.2)
    model = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, η)
    sol = solve(model)
    r = full_res(model, sol)
    g = real_gdp(sol)
    @printf("%6.2f  %.2e  %+.6f%%  %.6f  %.6f\n",
        η, max(maximum(abs.(r.zp)), maximum(abs.(r.mc)), abs(r.labor), abs(r.num)),
        100*(g/base_gdp - 1), sol.wages[1], r.debt)
end
println("DONE")
