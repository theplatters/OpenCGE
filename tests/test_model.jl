using BeyondHulten
using DataFrames
using Test

@testset "Baqaee--Farhi input-share convention" begin
	# Rows are suppliers and columns are users, as in the Destatis table.
	intermediate_use = [8.0 2.0; 3.0 7.0]
	Ω = BeyondHulten.conditional_input_shares(intermediate_use)

	@test Ω ≈ [8 / 11 3 / 11; 2 / 9 7 / 9]
	@test vec(sum(Ω, dims = 2)) ≈ ones(2)
	@test_throws ArgumentError BeyondHulten.conditional_input_shares([1.0 0.0; 2.0 0.0])
end

@testset "Impulse sector selection" begin
	io = DataFrame("Letzte Verwendung von Gütern zusammen" => fill(100.0, 71))
	data = (grossy = ones(71), io = io)
	impulses = DataFrame(zeros(2, 73), :auto) # year, 71 sectors, wages

	shocks = BeyondHulten.impulse_shock(data, impulses)
	@test length(shocks.demand_shock) == 71
	@test length(shocks.demand_shock_raw) == 71
	@test_throws DimensionMismatch BeyondHulten.impulse_shock(data, DataFrame(zeros(2, 72), :auto))
end
