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

@testset "Törnqvist real-GDP quantity index" begin
	base_prices = [1.0, 1.0]
	base_quantities = [1.0, 1.0]
	prices = [1.0, 2.0]
	quantities = [2.0, 1.0]

	index = tornqvist_quantity_index(prices, quantities, base_prices, base_quantities)
	@test index ≈ sqrt(2)
	@test tornqvist_quantity_index(3 .* prices, quantities, base_prices, base_quantities) ≈ index
	@test tornqvist_quantity_index(prices, quantities, base_prices, base_quantities) ≈
		tornqvist_quantity_index(base_prices, base_quantities, prices, quantities)^-1
	@test tornqvist_quantity_index([1.0, 2.0], [2.0, 0.0], [1.0, 1.0], [1.0, 0.0]) ≈ 2.0
	@test_throws DimensionMismatch tornqvist_quantity_index([1.0], [1.0, 2.0], [1.0], [1.0])
	@test_throws DomainError tornqvist_quantity_index([1.0, 1.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0])
end
