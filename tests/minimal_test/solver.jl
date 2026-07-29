# ── Newton solver with finite-difference Jacobian ──

function numerical_jacobian(f, x, params; eps=1e-6)
    n = length(x)
    f0 = f(x, params)
    m = length(f0)
    J = zeros(m, n)
    for j in 1:n
        x_plus = copy(x)
        x_plus[j] += eps
        f_plus = f(x_plus, params)
        J[:, j] = (f_plus - f0) / eps
    end
    return J
end

function solve_newton(f, x0, params; max_iter=500, tol=1e-10, verbose=false)
    x = copy(x0)
    n = length(x)

    for iter in 1:max_iter
        F = f(x, params)
        norm_F = norm(F)

        if verbose && (iter % 20 == 0 || iter == 1)
            @printf("    iter %d: ||F|| = %.2e\n", iter, norm_F)
        end

        if norm_F < tol
            return x, true, iter
        end

        J = numerical_jacobian(f, x, params)

        try
            dx = -(J \ F)
            alpha = min(1.0, 1.0 / max(1.0, norm(dx) / norm(x)))
            x_new = x + alpha * dx
            F_new = f(x_new, params)
            if norm(F_new) < norm_F
                x = x_new
            else
                x = x + 0.01 * alpha * dx
            end
        catch e
            x .+= 1e-4 * randn(n)
        end
    end

    return x, false, max_iter
end