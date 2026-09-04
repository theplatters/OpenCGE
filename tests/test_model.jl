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
    data = tiny_fixture()
    impulses = DataFrame(zeros(2, 4), :auto) # year, two sectors, wages

    shocks = BeyondHulten.impulse_shock(data, impulses)
    @test length(shocks.demand_shock) == length(data.grossy)
    @test length(shocks.demand_shock_raw) == length(data.grossy)
    @test_throws DimensionMismatch BeyondHulten.impulse_shock(data, DataFrame(zeros(2, 3), :auto))
end

@testset "shock construction and validation" begin
    data = tiny_fixture()
    for (constructor, field) in ((standard_shock, :demand_shock),
                               (standard_tech_shock, :supply_shock))
        shock = constructor(data, "b")
        @test length(getproperty(shock, field)) == 2
        @test getproperty(shock, field)[2] != 1
        @test getproperty(shock, field)[1] == 1
        @test_throws ArgumentError constructor(data, "missing")
    end
    shock = autonomous_shock(data; sector="b")
    @test length(shock.autonomous_demand) == 2
    @test shock.autonomous_demand == [0., 1.8097957577943152]
    @test_throws ArgumentError autonomous_shock(data; sector="missing")
    @test_throws DimensionMismatch Shocks(ones(2), ones(3), zeros(2))
    @test_throws DimensionMismatch Data(data.io, ones(2, 2), [.5], [.5, .5], [1., 1.],
        [.5, .5], [.5, .5], [1., 1.], [.5, .5])
    @test_throws DimensionMismatch Data(data.io, ones(1, 2), [.5, .5], [.5, .5], [1., 1.],
        [.5, .5], [.5, .5], [1., 1.], [.5, .5])
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
