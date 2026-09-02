# Diagnostic: run the package variance decomposition and descriptive η report.
include("bootstrap.jl")


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

println("\n=== eta_sweep_diagnostics (reduced grid) ===")
diagnostics = eta_sweep_diagnostics(data, shocks; θ=0.5, ϵ=0.5, σ=0.9)
summary_table(diagnostics.decomposition)
@printf("Maximum sectoral variation: %.6f\n", diagnostics.max_variation)
@printf("η first-order share: %.6f\n", diagnostics.eta_share)
