using ModelingToolkit
using BeyondHulten
using LinearAlgebra
using OrdinaryDiffEq
using Plots
using ModelingToolkit: t_nounits as t, D_nounits as D

function shock(t)
  1 + 0.0 * sin(2 * π * t)
end
@register_symbolic shock(t)
function create_economic_model(; sectors=2)
  # Define parameters
  @parameters begin
    σ = 1.01
    θ = 0.1
    ϵ = 0.1
    γ = 1.5
    ρ = -1.2
    δ = 0.1
    β = 1.001
    τ = 1.001
    Ω[1:2, 1:2] = [0.7 0.3; 0.5 0.5]
    A[1:2] = [1.0, 1.0]
    labor_share[1:2] = [0.2, 0.8]
  end

  # Collect all parameters
  ps = [σ, θ, ϵ, γ, ρ, δ, β, Ω, A, labor_share, τ]

  # Define variables with appropriate dimensions
  @variables begin
    output(t)[1:sectors] = ones(2)
    wages(t)[1:sectors] = ones(sectors)
    profits(t)[1:sectors]
    prices(t)[1:sectors] = ones(sectors)
    labor(t)[1:sectors]
    cost(t)[1:sectors]
    total_labor(t) = 1.0
    employment_ratio(t)
    X(t)[1:sectors, 1:sectors]
    consumption(t)[1:sectors]
    labor_productivity(t) = 1.0
    productivity(t)[1:sectors] = [1.0, 1.0]
  end
  # Create kquations
  eqs = Equation[]

  # Cost function calculations
  # We'll have to operate elementwise with ModelingToolkit
  for i in 1:sectors
    push!(eqs, productivity[i] ~ A[i] * shock(t))
    # Calculate cost for each sector
    labor_share_prod_i = labor_share[i] / labor_productivity
    G_i = sum(Ω[i, j] * prices[j]^(1 - θ) for j in 1:sectors)
    F_i = productivity[i]^(ϵ - 1) * (labor_share_prod_i) * wages[i]^(1 - ϵ) +
          (1 - labor_share_prod_i * G_i^((1 - ϵ) / (1 - θ)))
    costv_i = F_i^(1 / (1 - ϵ))

    # Calculate dcostdw for each sector
    dcostdw_i = productivity[i]^(ϵ - 1) * labor_share_prod_i * wages[i]^(-ϵ) * F_i^(ϵ / (1 - ϵ))

    # Calculate dcostdp for each sector and each price
    dcostdp_i = [productivity[i]^(ϵ - 1) * (1 - labor_share_prod_i) * Ω[i, j] *
                 prices[j]^(-θ) * G_i^(((1 - ϵ) / (1 - θ)) - 1) * F_i^(ϵ / (1 - ϵ))
                 for j in 1:sectors]

    # Dynamic equations
    push!(eqs, D(output[i]) ~ output[i] * profits[i] / (prices[i] * σ))
    push!(eqs, D(wages[i]) ~ (-γ + ρ * employment_ratio) + prices[i] * wages[i])
    push!(eqs, profits[i] ~ prices[i] * output[i] -
                            sum(X[i, j] * prices[j] for j in 1:sectors) -
                            wages[i] * labor[i])

    # Algebraic constraints
    push!(eqs, 0 ~ prices[i] - cost[i] * exp(δ * D(output[i])) - profits[i] / output[i])
    push!(eqs, 0 ~ cost[i] - costv_i)
    # Calculate term for labor constraint
    dcostdw_output_sum = sum(dcostdw_i * output[i] for i in 1:sectors)
    push!(eqs, 0 ~ labor[i] - min(total_labor / dcostdw_output_sum, 1) * dcostdw_i * output[i])
    push!(eqs, 0 ~ consumption[i] - output[i] - sum(X[i, j] for j in 1:sectors))
    # X matrix constraints
    for j in 1:sectors
      push!(eqs, 0 ~ X[i, j] - dcostdp_i[j] * output[i])
    end
  end

  # Total labor and employment equations
  push!(eqs, D(total_labor) ~ β * total_labor)
  push!(eqs, 0 ~ employment_ratio - sum(labor[i] for i in 1:sectors) / total_labor)

  push!(eqs, D(labor_productivity) ~ τ * labor_productivity)
  # Create an ODESystem with these equations
  @named sys = ODESystem(eqs, t,
    [output, wages, profits, prices, labor, cost,
      X, total_labor, employment_ratio],
    ps, tspan=(0.0, 1.0))

  return sys
end

# Create the model with 2 sectors
model = create_economic_model(sectors=2)

# Simplify the system for numerical solving
sys_simplified = structural_simplify(model, simplify=true, fully_determined=true)

# Set up the time span and solve
tspan = (0.0, 10.0)
prob = ODEProblem(sys_simplified, [], tspan, [])
sol = OrdinaryDiffEq.solve(prob, RadauIIA9())
plot(sol, size=(1980, 1080))
[x[1] for x in sol[sys_simplified.profits]] |> plot
[sum(x) for x in sol[sys_simplified.X]] |> plot
sol[sys_simplified.profits]
sol[sys_simplified.employment_ratio] |> plot
sol[sys_simplified.wages]
sol[sys_simplified.cost]
sol[sys_simplified.prices]

