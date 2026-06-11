using Test

include("main.jl")
include("test_model.jl")

@testset "All Tests" begin
    include("test_model.jl")
end