# Comprehensive bridge gauge: test consumption-Törnqvist under multiple shocks.
# Uses proper no-shock baseline for all Δ% calculations.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
N = 71
SEC = 35  # construction

# Proper no-shock baseline
ref = mobile_labor_model(data, Shocks(ones(N), ones(N), zeros(N)), 0.5, 2.0, 0.9, 0.0)
ref_sol = solve(ref)
ref_gdp = real_gdp(ref_sol)
@printf("No-shock baseline real_gdp = %.10f\n\n", ref_gdp)

function bridge(model_fn, label)
    etas = [0.0, 0.25, 0.5, 0.75, 1.0]
    gdps = Float64[]
    resids = Float64[]
    for η in etas
        m = model_fn(η)
        s = solve(m)
        push!(gdps, real_gdp(s))
        push!(resids, maximum(abs, equilibrium_residuals(s)))
    end
    d0 = 100*(gdps[1]/ref_gdp - 1)
    d1 = 100*(gdps[end]/ref_gdp - 1)
    br = d1 - d0
    @printf("%-54s  η=0: %+.6f%%  η=1: %+.6f%%  bridge: %s%+.6f%%  maxRes: %.2e\n",
        label, d0, d1, br >= 0 ? "+" : "", br, maximum(resids))
    return (d0, d1, br)
end

println("═"^78)
println("BRIDGE: immobile (η=0) vs mobile (η=1) — consumption Törnqvist (B&F metric)")
println("═"^78)

# 1. Supply shock: +30% to construction (ε=2, substitution)
bridge(
    η -> mobile_labor_model(data, Shocks(begin ss=ones(N); ss[SEC]*=1.30; ss end, ones(N), zeros(N)), 0.5, 2.0, 0.9, η),
    "A  Supply +30% (ε=2.0)"
)

# 2. Supply shock: +100% to construction (ε=2)
bridge(
    η -> mobile_labor_model(data, Shocks(begin ss=ones(N); ss[SEC]*=2.0; ss end, ones(N), zeros(N)), 0.5, 2.0, 0.9, η),
    "B  Supply +100% (ε=2.0)"
)

# 3. Supply shock: +30% (ε=0.5, complementarity)
bridge(
    η -> mobile_labor_model(data, Shocks(begin ss=ones(N); ss[SEC]*=1.30; ss end, ones(N), zeros(N)), 0.5, 0.5, 0.9, η),
    "C  Supply +30% (ε=0.5)"
)

# 4. Demand shock: consumption shift +30% to construction (ε=0.5)
bridge(
    η -> mobile_labor_model(data, Shocks(ones(N), begin ds=ones(N); ds[SEC]*=1.30; ds end, zeros(N)), 0.5, 0.5, 0.9, η),
    "D  Demand-shift +30% (ε=0.5)"
)

# 5. Demand shock: consumption shift +30% to construction (ε=2.0)
bridge(
    η -> mobile_labor_model(data, Shocks(ones(N), begin ds=ones(N); ds[SEC]*=1.30; ds end, zeros(N)), 0.5, 2.0, 0.9, η),
    "E  Demand-shift +30% (ε=2.0)"
)

# 6. Autonomous demand (construction, mult=0.2, ε=0.5) — small, stable
bridge(
    η -> mobile_labor_model(data, autonomous_shock(data; autonomous_mult=0.2), 0.5, 0.5, 0.9, η),
    "F  Autonomous demand 0.2× (ε=0.5)"
)

# 7. Change sigma — lower σ = less consumption substitution
bridge(
    η -> mobile_labor_model(data, Shocks(begin ss=ones(N); ss[SEC]*=1.30; ss end, ones(N), zeros(N)), 0.5, 2.0, 0.5, η),
    "G  Supply +30% (ε=2.0, σ=0.5)"
)

# 8. Change sigma — very low σ
bridge(
    η -> mobile_labor_model(data, Shocks(begin ss=ones(N); ss[SEC]*=1.30; ss end, ones(N), zeros(N)), 0.5, 2.0, 0.1, η),
    "H  Supply +30% (ε=2.0, σ=0.1)"
)

println("─"^78)
println("Interpretation: the reallocation bridge is a second-order Harberger effect.")
println("B&F (2019) find it quantitatively tiny in competitive one-factor models;")
println("the original '+1.5%' was an artifact of the buggy non-equilibrium.")
println("DONE")
