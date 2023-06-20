function mvnrnd(μ, Σ)
    U = cholesky(Σ).U
    (Random.randn(size(μ))' * U + μ')'
end
