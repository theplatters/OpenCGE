# Full descriptive η sweep and Sobol decomposition through the package API.
include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
diagnostics = eta_sweep_diagnostics(data, standard_shock(data))

println("η        real GDP      wage")
for (η, sol) in zip(diagnostics.sweep.η_values, diagnostics.sweep.solutions)
    @printf("%-8.2f %.6f    %.6f\n", η, real_gdp(sol), sol.wages_raw[1])
end

summary_table(diagnostics.decomposition)
@printf("Maximum sectoral variation: %.6f\n", diagnostics.max_variation)
@printf("η first-order share: %.6f\n", diagnostics.eta_share)
