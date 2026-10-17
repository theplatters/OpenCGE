using OpenCGE, Test, LinearAlgebra

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the final-demand bundle and the three welfare indices
# (2026-10-16): `final_demand_bundle`, `ces_quantity_index`,
# `real_consumption_ev`, `real_final_demand_ces`, `welfare_leisure` in
# src/core/diagnostics.jl. Specification: revision/proposals/45-welfare-indices.md.
#
# The load-bearing canary here is C10: a CES quantity index must NOT equal the
# value ratio for a basket that does not grow proportionally. The earlier
# specification asserted the opposite and got the sign of a reported cell wrong.

@testset "C10 lopsided basket: CES index is not the value ratio" begin
	# bread and cement: 100/10 -> 95/30. Prices unchanged at 1.
	q0 = [100.0, 10.0]
	q1 = [95.0, 30.0]
	vr = sum(q1) / sum(q0)
	@test vr ≈ 1.136363636 rtol = 1e-9
	# CES at the paper's sigma: complements, so the lopsided basket is penalised
	@test ces_quantity_index(q1, q0, 0.9) ≈ 1.048520246 rtol = 1e-8
	@test ces_quantity_index(q1, q0, 0.5) ≈ 1.012924071 rtol = 1e-8
	@test ces_quantity_index(q1, q0, 2.0) ≈ 1.088957736 rtol = 1e-8
	# every sigma < 1 lands far below the value ratio; sigma = 2 still does
	@test ces_quantity_index(q1, q0, 0.9) < vr
	@test ces_quantity_index(q1, q0, 2.0) < vr
	# the Tornqvist is share-weighted and moves the other way --- which is why
	# "the indices agree" cannot be rescued by swapping the functional form
	w = q0 ./ sum(q0)
	@test exp(sum(((q1 ./ sum(q1)) .+ w) ./ 2 .* log.(q1 ./ q0))) ≈ 1.149082515 rtol = 1e-8
end

@testset "ces_quantity_index: exact cases and guards" begin
	q0 = [100.0, 10.0]
	# proportional growth: every index, at every sigma, equals the value ratio
	for s in (0.5, 0.9, 1.0, 2.0)
		@test ces_quantity_index(2.0 .* q0, q0, s) ≈ 2.0 rtol = 1e-12
		@test ces_quantity_index(0.5 .* q0, q0, s) ≈ 0.5 rtol = 1e-12
	end
	# sigma = 1 is the Cobb-Douglas (geometric) index
	q1 = [95.0, 30.0]
	@test ces_quantity_index(q1, q0, 1.0) ≈ exp(sum((q0 ./ sum(q0)) .* log.(q1 ./ q0))) rtol = 1e-12
	# guards
	@test_throws DimensionMismatch ces_quantity_index([1.0], q0, 0.9)
	@test_throws DomainError ces_quantity_index([-1.0, 2.0], q0, 0.9)
	@test_throws DomainError ces_quantity_index(q0, zeros(2), 0.9)
	# a component with a zero baseline is dropped, not divided by zero
	@test isfinite(ces_quantity_index([95.0, 30.0], [100.0, 0.0], 0.9))
end

@testset "final-demand bundle and indices on the open-economy fixture" begin
	data = tiny_fixture()
	N = length(data.factor_share)
	sh0 = Shocks(ones(N), ones(N), zeros(N))
	g = [0.007, 0.003]
	base_m = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0)
	base = solve(base_m)
	m = mobile_labor_model(data, sh0, 0.5, 0.5, 0.9, 1.0; financing = TaxFinanced(g))
	sol = solve(m; init = [ones(N); data.λ; 1.0])

	b = final_demand_bundle(m, sol.prices_raw, sol.quantities, OpenCGE._bundle_wage(sol);
		external_transfer = sol.external_transfer)
	@test length(b.c) == N && length(b.g) == N && length(b.f) == N
	@test b.f ≈ b.c .+ b.g                      # lambda = 1 default
	b0 = final_demand_bundle(m, sol.prices_raw, sol.quantities, OpenCGE._bundle_wage(sol);
		lambda = 0.0, external_transfer = sol.external_transfer)
	@test b0.f ≈ b0.c

	# baseline normalisation
	@test real_final_demand_ces(base, base) ≈ 1.0 rtol = 1e-12
	@test real_consumption_ev(base, base) ≈ 1.0 rtol = 1e-12
	@test welfare_leisure(base, base).welfare ≈ 0.0 atol = 1e-12

	# C9: lambda = 0 returns the pure consumption measure
	@test real_final_demand_ces(sol, base; lambda = 0.0, sigma = 0.9) ≈
		real_consumption_ev(sol, base) rtol = 1e-6

	# the programme moves the public bundle, so lambda changes the answer
	q1 = real_final_demand_ces(sol, base; lambda = 1.0)
	q0 = real_final_demand_ces(sol, base; lambda = 0.0)
	@test q1 > q0

	# C7/C8 scope: the value-ratio variant, not the CES index, is what satisfies
	# the wage-pinned accounting identities; both are returned so a caller can
	# tell which one is being quoted.
	wl = welfare_leisure(sol, base; lambda = 1.0)
	@test haskey(wl, :value_ratio) && haskey(wl, :q) && haskey(wl, :leisure_cost)
	@test wl.value_ratio != wl.q

	# the saving convention enters only through the scale of the private part
	lk = welfare_leisure(sol, base; saving = :leak)
	df = welfare_leisure(sol, base; saving = :deferred)
	@test isfinite(lk.welfare) && isfinite(df.welfare)
	# the compact synthetic fixture carries no saving (s = 0), so the two
	# conventions must coincide exactly here; a non-zero-saving case is
	# exercised paper-side by welfare-fd-ces.py and the v10 manifests
	@test data.saving_rate == 0.0
	@test df.welfare ≈ lk.welfare rtol = 1e-12
	# the mechanism itself: stretching the private part changes the CES index
	# whenever the basket is not proportional --- which is the general case
	@test ces_quantity_index(1.12 .* b.c .+ b.g, b.f, 0.9) != ces_quantity_index(b.c .+ b.g, b.f, 0.9)

	# C11: the own-weights aggregator IS the bundle's value ratio, so the
	# leisure-augmented form on it equals (value_ratio - 1) - leisure_cost
	# identically. This is the measure the bundle column and W_c carry from v6
	# on; the fixed-weight CES stays available as the stanced alternative.
	ow = welfare_leisure(sol, base; lambda = 1.0, aggregator = :ownweights)
	@test ow.q ≈ ow.value_ratio rtol = 1e-12
	@test ow.welfare ≈ (ow.value_ratio - 1) - ow.leisure_cost rtol = 1e-12
	@test welfare_leisure(base, base; aggregator = :ownweights).welfare ≈ 0.0 atol = 1e-12
	# the fixture's basket does not grow proportionally, so the two aggregators
	# must give different answers -- if they agreed, v6 would be a no-op
	@test ow.welfare != wl.welfare
	@test ow.welfare == welfare_leisure(sol, base; aggregator = :ces).welfare ||
		ow.welfare != wl.welfare

	# errors
	@test_throws ArgumentError welfare_leisure(sol, base; saving = :nope)
	@test_throws ArgumentError welfare_leisure(sol, base; aggregator = :nope)
	@test_throws ArgumentError final_demand_bundle(Model(data, sh0, CES()),
		sol.prices_raw, sol.quantities, 1.0)
end
