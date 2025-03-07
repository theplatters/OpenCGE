using BeyondHulten
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using LinearAlgebra

data = Data("I-O_DE2019_formatiert.csv"; read_unemployment=false)
good_elasticity = (0.9 - 1) / 0.9
good_labor_elasticity = (0.9 - 1) / 0.9
ces(x, t, l, A) = A * (((1 .- data.factor_share) .* diag(data.Ω * x .^ good_elasticity)) .^  (good_labor_elasticity/good_elasticity)  + data.factor_share .* l .^ good_labor_elasticity) .^ (1 / good_labor_elasticity)

ces(ones(71, 71), 1, ones(71), 1)
@mtkmodel GOODWIN_CES begin
  @parameters begin
    α = 0.5    # Productivity growth rate
    β = 0.4    # Labor productivity parameter
    γ = 1.6    # Capital depreciation rate
    ρ = 3.4    # Output-capital ratio
    σ = 1.3    # Original model parameter
    
    # CES-specific parameters
    η = 0.7         # Capital factor share (1-η = labor share)
    Ω = 1.2         # Capital efficiency parameter
    A = 1.0         # Technology level
    σ_ces = 0.8     # Elasticity of substitution (σ = 1/(1-ρ_ces))
  end
  
  @variables begin
    v(t) = 1.0      # Wage share
    u(t) = 1.0      # Employment rate
    K(t) = 2.0      # Capital stock (new variable)
  end

  @equations begin
    # CES Production Function (simplified form)
    ρ_ces = 1 - 1/σ_ces  # Derived substitution parameter
    Y ~ A * ((η * (Ω*K)^ρ_ces + (1-η) * u^ρ_ces))^(1/ρ_ces)

    # Modified Goodwin equations incorporating CES production
    D(v) ~ v * ( (Y - β * u)/Y - α - γ )  # Wage share dynamics
    D(u) ~ u * (ρ * (Y/K) - α - β)        # Employment dynamics
    D(K) ~ Y - γ*K                        # Capital accumulation
  end
end

import OrdinaryDiffEq as OD
@mtkbuild fol = GOODWIN_CES()
prob = OD.ODEProblem(fol, [], (0.0, 100.0), [])
sol = OD.solve(prob)

using Plots
plot(sol)
