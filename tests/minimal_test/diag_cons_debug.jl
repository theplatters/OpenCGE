# Debug: check what's happening with consumption and Törnqvist at η=0 under supply shock
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
N = 71
SEC = 35

# Baseline (no shock)
base = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 2.0, 0.9, 0.0)
bs = solve(base)
@printf("Baseline real_gdp = %.10f\n", real_gdp(bs))
@printf("Baseline sum(consumption) = %.10f\n", sum(bs.consumption))
@printf("Baseline wage = %.10f\n", bs.wages[1])
@printf("Baseline p[35] = %.10f\n", bs.prices[35])
println()

# η=0 under +30% supply shock
ss = ones(N); ss[SEC] *= 1.30
m0 = mobile_labor_model(data, Shocks(ss, ones(N), zeros(N)), 0.5, 2.0, 0.9, 0.0)
s0 = solve(m0)
@printf("η=0 shocked real_gdp = %.10f\n", real_gdp(s0))
@printf("η=0 Δreal_gdp = %+.10f%%\n", 100*(real_gdp(s0)/real_gdp(bs) - 1))
@printf("η=0 sum(consumption) = %.10f\n", sum(s0.consumption))
@printf("η=0 wage = %.10f\n", s0.wages[1])
@printf("η=0 p[35] = %.10f\n", s0.prices[35])

# Manually compute Törnqvist using the function
rgdp = tornqvist_quantity_index(
    s0.prices, s0.consumption,
    ones(N), data.consumption_share .* sum(data.labor_share)
)
println("\nManual Törnqvist = ", rgdp)

# Check: what are the consumption vectors?
println("\nTop 5 consumption (baseline):")
for i in sortperm(bs.consumption, rev=true)[1:5]
    @printf("  [%2d] shocked=%.6f  base_c=%.6f  share=%.6f  p=%.6f\n",
        i, s0.consumption[i], bs.consumption[i], data.consumption_share[i], s0.prices[i])
end
println("consumption_share[35] = ", data.consumption_share[SEC])
@printf("consumption[35] base=%.6f shocked=%.6f\n", bs.consumption[SEC], s0.consumption[SEC])

# What is the Tornqvist log growth contribution for sector 35?
base_cons = data.consumption_share .* sum(data.labor_share)
exp_base = ones(N) .* base_cons
exp_cur = s0.prices .* s0.consumption
total_base = sum(exp_base)
total_cur = sum(exp_cur)
sb = exp_base ./ total_base
sc = exp_cur ./ total_cur
active = (sb .+ sc) .> 0
@printf("\nSector 35: share_base=%.8f share_cur=%.8f log(q_cur/q_base)=%.8f  contrib=%.8f\n",
    sb[SEC], sc[SEC], log(s0.consumption[SEC]/base_cons[SEC]),
    0.5*(sb[SEC]+sc[SEC])*log(s0.consumption[SEC]/base_cons[SEC]))
println("\nDONE")