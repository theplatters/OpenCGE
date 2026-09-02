# Diagnostic: does the BASE CES model (sector-specific labor, full_labor_slack)
# also invert under a positive demand shock? Determines if the bug is
# specific to mobile_labor or systemic.
include("bootstrap.jl")


data = Data("I-O_DE2019_formatiert.csv")
N = length(data.factor_share)

# Baseline (no shock)
base = Shocks(ones(N), ones(N), zeros(N))
m0 = Model(data, base, CES(CESElasticities(0.5,0.5,0.9), full_labor_slack, false))
s0 = solve(m0)
@printf("BASELINE  real_gdp=%.6f  nominal_gdp=%.6f  mean_wage=%.6f\n", real_gdp(s0), nominal_gdp(s0), mean(s0.wages))

# Standard +81% demand shock in construction
shocks = standard_shock(data)
m1 = Model(data, shocks, CES(CESElasticities(0.5,0.5,0.9), full_labor_slack, false))
s1 = solve(m1)
@printf("SHOCKED   real_gdp=%.6f  nominal_gdp=%.6f  mean_wage=%.6f\n", real_gdp(s1), nominal_gdp(s1), mean(s1.wages))
@printf("  -> real GDP change vs baseline: %+.4f%%\n", 100*(real_gdp(s1)-real_gdp(s0)))
@printf("  -> wage change vs baseline:     %+.4f%%\n", 100*(mean(s1.wages)-mean(s0.wages))/mean(s0.wages))
