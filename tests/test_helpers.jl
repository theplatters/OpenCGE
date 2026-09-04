using DataFrames
using LinearAlgebra

"""Small, self-contained two-sector calibration used throughout the tests."""
function tiny_fixture()
    io = DataFrame("Sektoren" => ["a", "b"],
        "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0])
    Data(io, Matrix{Float64}(I, 2, 2), [.5, .5], [.5, .5], [1., 1.],
        [.5, .5], [.5, .5], [1., 1.], [0.5, 0.5])
end

max_equilibrium_residual(sol) = maximum(abs, equilibrium_residuals(sol))
