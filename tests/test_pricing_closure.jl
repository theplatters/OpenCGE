using OpenCGE, DataFrames, LinearAlgebra, Test

isdefined(Main, :tiny_fixture) || include(joinpath(@__DIR__, "test_helpers.jl"))

# Contract tests for the pricing door (ADR-0029, Foundation III): the
# capacity / Verdoorn route (CapacityPricing) and the markup route
# (MarkupPricing). Synthetic fixtures only. The open fixture carries saving,
# an import margin and injections (s > 0), so the fixed-wage η = 1 system is
# determinate; `split_bill` hand-splits its purchaser-price bill
# (A_bill + M_int + T_int = (1 − fs)·λ) so the bill-factor sites of the canary
# and `gdp_components` are actually exercised (on the compact constructor
# M_int = T_int = 0 and the identity stays green whatever those sites compute).

"""Open v3-style 3-sector fixture (saving, import margin, gov/inv/exports)."""
function pricing_fixture()
	io = DataFrame("Sektoren" => ["a", "b", "c"],
		"Letzte Verwendung von Gütern zusammen" => [1.0, 1.0, 1.0])
	Ω = [0.6 0.3 0.1; 0.2 0.5 0.3; 0.25 0.25 0.5]
	fs = [0.6, 0.5, 0.4]
	grossy = [2.0, 1.5, 1.0]
	va = fs .* grossy
	gdp = sum(va)
	λ = grossy ./ gdp
	m = 0.2
	gov = [0.03, 0.015, 0.005]
	inv = [0.03, 0.075, 0.045]
	expo = [0.075, 0.045, 0.03]
	c0 = (λ .- Ω' * ((1 .- fs) .* λ) .- (1 - m) .* (gov .+ inv) .- expo) ./ (1 - m)
	s = 1 - sum(c0) / (1 - sum(gov))
	ω = c0 ./ sum(c0)
	Data(io, Ω, Ω, ω, fs, λ, λ .* fs, ω, grossy, va, grossy, DataFrame(),
		zeros(3), zeros(3), zeros(3), gov, c0, fill(m, 3), inv, expo, s,
		(1 .- fs) .* λ, zeros(3), zeros(3), gdp, gdp, gdp)
end

"""`d` with its purchaser-price bill split 70/20/10 into A_bill/M_int/T_int."""
function split_bill(d::Data)
	vals = Any[getfield(d, f) for f in fieldnames(Data)]
	bill = (1 .- d.factor_share) .* d.λ
	vals[findfirst(==(:A_bill), fieldnames(Data))] = 0.7 .* bill
	vals[findfirst(==(:M_int), fieldnames(Data))] = 0.2 .* bill
	vals[findfirst(==(:T_int), fieldnames(Data))] = 0.1 .* bill
	Data(vals...)
end

const _PR_G = [0.02, 0.0, 0.0]
_pr_demand() = Shocks(ones(3), [1.2, 1.0, 1.0], zeros(3))
_pr_supply() = Shocks([1.1, 1.0, 0.95], [1.2, 1.0, 1.0], zeros(3))
_pr_vec(m, s) = OpenCGE._solution_vector(m, s)

# (closure, η, financing) cells shared by the nesting and account testsets.
const _PR_CELLS = [
	(:fixed, 1.0, TaxFinanced(_PR_G)), (:fixed, 1.0, ExternalDebt(_PR_G)),
	(:mobile, 1.0, TaxFinanced(_PR_G)), (:mobile, 1.0, ExternalDebt(_PR_G)),
	(:mobile, 0.0, TaxFinanced(_PR_G)), (:mobile, 0.0, ExternalDebt(_PR_G)),
]

"""Same solution, to the last bit, in every reported field."""
function bit_identical(m0, s0, m1, s1)
	_pr_vec(m0, s0) == _pr_vec(m1, s1) && s0.consumption == s1.consumption &&
		s0.real_gdp == s1.real_gdp && s0.nominal_gdp == s1.nominal_gdp &&
		equilibrium_residuals(s0) == equilibrium_residuals(s1) &&
		gdp_components(m0, s0) == gdp_components(m1, s1)
end

@testset "pricing: types and NoPricing defaults" begin
	d = pricing_fixture()
	sh = _pr_demand()
	@test Model(d, sh, MobileLaborCES(MobileLaborCESElasticities(0.5, 0.5, 0.9, 1.0), d)).pricing isa NoPricing
	@test Model(d, sh, MobileLaborCES(MobileLaborCESElasticities(0.5, 0.5, 0.9, 1.0), d), NoFinancing()).pricing isa NoPricing
	@test mobile_labor_model(d, sh, 0.5, 0.5, 0.9, 1.0).pricing isa NoPricing
	@test mobile_labor_model(d, sh, 0.5, 0.5, 0.9, 1.0; pricing = CapacityPricing(0.2)).pricing == CapacityPricing(0.2)
	@test_throws ArgumentError CapacityPricing(Inf)
	@test_throws ArgumentError MarkupPricing(NaN)
	@test_throws ArgumentError MarkupPricing(0.5; s_pi = 1.5)
	@test MarkupPricing(0.5).s_pi == 1.0
	# NoPricing hooks return their inputs (the same objects, not copies).
	A, y, λ, p = [1.1, 1.0, 0.9], [1.0, 2.0, 3.0], [1.0, 1.0, 1.0], [1.0, 1.1, 0.9]
	@test pricing_technology(NoPricing(), A, y, λ) === A
	@test pricing_markup(NoPricing(), y, λ) === nothing
	@test pricing_foc_price(NoPricing(), p, nothing) === p
	@test pricing_profit(NoPricing(), p, nothing, y) === nothing
end

# ── Route 1: capacity / Verdoorn ─────────────────────────────────────────────

@testset "pricing: δ = 0 nests the baseline bit for bit" begin
	d = pricing_fixture()
	for data in (d, split_bill(d)), sh in (_pr_demand(), _pr_supply()),
			(cl, η, fin) in _PR_CELLS
		m0 = mobile_labor_model(data, sh, 0.5, 0.5, 0.9, η; closure = cl, financing = fin)
		m1 = mobile_labor_model(data, sh, 0.5, 0.5, 0.9, η; closure = cl, financing = fin,
			pricing = CapacityPricing(0.0))
		@test bit_identical(m0, solve(m0), m1, solve(m1))
	end
end

@testset "pricing: demand-only shock moves prices only through the door" begin
	# Lemma 1 of the paper: under A ≡ 1 the zero-profit block pins p = 1 in
	# every single-wage closure. The capacity door breaks it.
	d = pricing_fixture()
	for (cl, η, fin) in _PR_CELLS[1:4]
		m0 = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, η; closure = cl, financing = fin)
		@test maximum(abs, solve(m0).prices .- 1) ≤ 1e-10
		for δ in (0.5, -0.5)
			m = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, η; closure = cl,
				financing = fin, pricing = CapacityPricing(δ))
			s = solve(m)
			@test maximum(abs, equilibrium_residuals(s)) ≤ 1e-9
			@test maximum(abs, s.prices .- 1) > 1e-2   # measured 0.023–0.052
		end
	end
end

@testset "pricing: Cobb–Douglas price identity (capacity)" begin
	# θ = ϵ = 1, w = 1: log p = (I − B)⁻¹ (log μ − log A_eff), B = diag(1−fs)·Ω_raw,
	# with μ ≡ 1 on the capacity route and A_eff = A (y/λ)^(−δ).
	d = pricing_fixture()
	B = Diagonal(1 .- d.factor_share) * d.Ω_raw
	for δ in (0.5, -0.5)
		m = mobile_labor_model(d, _pr_supply(), 1.0, 1.0, 0.9, 1.0; closure = :fixed,
			financing = TaxFinanced(_PR_G), pricing = CapacityPricing(δ))
		s = solve(m)
		A_eff = m.shocks.supply_shock .* (s.quantities ./ d.λ) .^ (-δ)
		@test log.(s.prices_raw) ≈ (I - B) \ (-log.(A_eff)) atol = 1e-9
	end
end

@testset "pricing: external account closes on the hand-split bill (capacity)" begin
	d = split_bill(pricing_fixture())
	@test d.A_bill .+ d.M_int .+ d.T_int ≈ (1 .- d.factor_share) .* d.λ atol = 1e-15
	for (cl, η, fin) in _PR_CELLS, δ in (0.5, -0.5)
		m = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = CapacityPricing(δ))
		s = solve(m)
		X = _pr_vec(m, s)
		can = external_balance_canary(m, X)
		@test abs(can.diff) ≤ 1e-9
		@test gdp_components(m, s).wedge ≈ -can.diff atol = 1e-12
		# The fixture exercises the bill-factor site: valuing the leaks with
		# the raw supply shock instead of A_eff opens the identity.
		ϵ = m.options.elasticities.ϵ
		ip = OpenCGE._intermediate_price(d.Ω_raw, s.prices_raw, m.options.elasticities.θ)
		k = pricing_bill_factor(s.prices_raw, m.shocks.supply_shock .* (s.quantities ./ d.λ) .^ (-δ), ip, ϵ)
		k_raw = pricing_bill_factor(s.prices_raw, m.shocks.supply_shock, ip, ϵ)
		@test abs(dot((k_raw .- k) .* ((d.M_int .+ d.T_int) ./ d.λ), s.quantities)) > 1e-4
	end
end

@testset "pricing: capacity guard and uniqueness" begin
	# The ADR-0017 guard runs on a NoPricing clone, so a closed fixture is
	# still rejected at the fixed-wage η = 1 system.
	data = tiny_fixture()
	@test_throws ArgumentError solve(mobile_labor_model(data, Shocks(ones(2), ones(2), zeros(2)),
		0.5, 0.5, 0.9, 1.0; closure = :fixed, pricing = CapacityPricing(0.3)))
	# Increasing returns (δ < 0) is positive feedback; the multi-start runs
	# agree on the open fixture (s > 0). Measured spread ≤ 1.7e-10, cond ≈ 20–29.
	d = pricing_fixture()
	for (θ, ϵ) in ((0.5, 0.5), (1.0, 1.0))
		m = mobile_labor_model(d, _pr_demand(), θ, ϵ, 0.9, 1.0; closure = :fixed,
			financing = TaxFinanced(_PR_G), pricing = CapacityPricing(-0.5))
		s = solve(m)
		adm = pricing_admissibility(m, s)
		@test adm.converged == length(adm.starts) == 4
		@test adm.spread ≤ 1e-8
		@test adm.rank == adm.n == 6
		@test 1 ≤ adm.cond < 1e6
	end
end

# ── Route 2, Stage A: the retained markup ────────────────────────────────────

@testset "pricing: κ = 0 nests the baseline bit for bit" begin
	d = pricing_fixture()
	for data in (d, split_bill(d)), sh in (_pr_demand(), _pr_supply()),
			(cl, η, fin) in _PR_CELLS
		m0 = mobile_labor_model(data, sh, 0.5, 0.5, 0.9, η; closure = cl, financing = fin)
		m1 = mobile_labor_model(data, sh, 0.5, 0.5, 0.9, η; closure = cl, financing = fin,
			pricing = MarkupPricing(0.0))
		s0, s1 = solve(m0), solve(m1)
		@test bit_identical(m0, s0, m1, s1)
		@test external_balance_canary(m0, _pr_vec(m0, s0)).diff ==
			external_balance_canary(m1, _pr_vec(m1, s1)).diff
	end
end

@testset "pricing: labour FOC at marginal cost (Shephard)" begin
	# At a κ > 0 solution, y·∂c/∂w = L_i with c the MARGINAL cost (p/μ). The
	# FOC evaluated at the purchaser price p inflates every L_i by μ^ϵ.
	d = pricing_fixture()
	m = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
		financing = TaxFinanced(_PR_G), pricing = MarkupPricing(0.5))
	s = solve(m)
	p, y = s.prices_raw, s.quantities
	(; ϵ, θ) = m.options.elasticities
	ip = OpenCGE._intermediate_price(d.Ω_raw, p, θ)
	c(w) = OpenCGE._ces_unit_cost(m.shocks.supply_shock, d.factor_share, w, ip, ϵ)
	h = 1e-6
	shephard = (c(1 + h) .- c(1 - h)) ./ (2h) .* y
	L = sectoral_labor_demand(p, y, 1.0, m)
	@test L ≈ shephard rtol = 1e-7
	@test p ≈ pricing_markup(m.pricing, y, d.λ) .* c(1.0) rtol = 1e-10
	m_p = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
		financing = TaxFinanced(_PR_G))   # same (p, y), FOC at p
	@test maximum(abs, sectoral_labor_demand(p, y, 1.0, m_p) .- shephard) > 1e-3
end

@testset "pricing: Stage A leaves the household budget unchanged" begin
	d = pricing_fixture()
	for (cl, η, fin) in _PR_CELLS[1:4]
		m = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(0.5))
		s = solve(m)
		w = s.wages_raw[1]
		L = sum(sectoral_labor_demand(s.prices_raw, s.quantities, w, m))
		τ = tau_rate(fin, m, s.prices_raw, w, L)
		@test household_expenditure(fin, m, w * L, s.prices_raw, L;
			external_transfer = s.external_transfer) ≈
			(1 - τ) * w * L + s.external_transfer atol = 1e-14
		@test OpenCGE._solution_blocks(m, s)[2].E ≈
			(1 - τ) * w * L + s.external_transfer atol = 1e-14
	end
end

@testset "pricing: external account closes with the retained rent booked" begin
	d = split_bill(pricing_fixture())
	for (cl, η, fin) in _PR_CELLS, κ in (0.25, 0.5)
		m = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(κ))
		s = solve(m)
		can = external_balance_canary(m, _pr_vec(m, s))
		Π = retained_profit(m, s)
		@test abs(Π) > 1e-3                           # the cell carries a rent
		@test can.retained_profit ≈ Π atol = 1e-14
		@test abs(can.diff) ≤ 1e-9                    # booked: the identity closes
		@test can.diff - can.retained_profit ≈ -Π atol = 1e-9   # unbooked: short by Π
		gc = gdp_components(m, s)
		@test gc.profit ≈ Π atol = 1e-14
		@test gc.income ≈ gc.wage_bill + Π atol = 1e-14
		@test gc.wedge ≈ -can.diff atol = 1e-12
	end
	m0 = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, 1.0)
	@test retained_profit(m0, solve(m0)) == 0.0
	@test external_balance_canary(m0, _pr_vec(m0, solve(m0))).retained_profit == 0.0
end

@testset "pricing: legacy-manna reading unchanged by the rent term" begin
	# ADR-0005 path: the gap reads the unbooked manna p·(A+G), with the rent
	# booked, exactly as without a markup.
	d = pricing_fixture()
	autonomous = [0.1, 0.0, 0.0]
	sh = Shocks(ones(3), ones(3); autonomous_demand = autonomous)
	A = autonomous .* d.consumption_share .* sum(d.labor_share)
	for pr in (NoPricing(), MarkupPricing(0.5))
		m = mobile_labor_model(d, sh, 0.5, 0.5, 0.9, 1.0; pricing = pr)
		s = solve(m)
		can = external_balance_canary(m, _pr_vec(m, s))
		@test dot(s.prices_raw, A) > 1e-3
		@test can.diff ≈ dot(s.prices_raw, A) atol = 1e-9
		pr isa MarkupPricing && @test abs(can.retained_profit) > 1e-5
	end
end

@testset "pricing: monotone in κ" begin
	# Demand-only shock on sector 1 (fixed wage, F2): the expanding sector's
	# price rises and the contracting sectors' fall as κ grows.
	d = pricing_fixture()
	sols = map(0.0:0.1:0.5) do κ
		solve(mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
			financing = TaxFinanced(_PR_G), pricing = MarkupPricing(κ)))
	end
	p1 = [s.prices[1] for s in sols]
	p2 = [s.prices[2] for s in sols]
	@test p1[1] ≈ 1.0 atol = 1e-12
	@test all(diff(p1) .> 0)
	@test all(diff(p2) .< 0)
	@test all(diff([s.real_gdp for s in sols]) .< 0)
end

@testset "pricing: Cobb–Douglas price identity (markup)" begin
	# The capacity identity with log μ ≠ 0; contains the uniform-markup
	# cascade p = μ̄^((I−B)⁻¹ 1) as a special case (not p = μ̄).
	d = pricing_fixture()
	B = Diagonal(1 .- d.factor_share) * d.Ω_raw
	for κ in (0.25, 0.5)
		m = mobile_labor_model(d, _pr_supply(), 1.0, 1.0, 0.9, 1.0; closure = :fixed,
			financing = TaxFinanced(_PR_G), pricing = MarkupPricing(κ))
		s = solve(m)
		μ = pricing_markup(m.pricing, s.quantities, d.λ)
		@test log.(s.prices_raw) ≈ (I - B) \ (log.(μ) .- log.(m.shocks.supply_shock)) atol = 1e-9
	end
end

@testset "pricing: non-positive markup equilibrium raises a named error" begin
	d = pricing_fixture()
	m = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
		financing = TaxFinanced(_PR_G), pricing = MarkupPricing(0.5))
	X = _pr_vec(m, solve(m))
	@test OpenCGE._assert_pricing_feasible(m, m.pricing, X) === nothing
	Xneg = copy(X); Xneg[4] = -0.01                     # y_1 < 0
	@test_throws PricingFeasibilityError OpenCGE._assert_pricing_feasible(m, m.pricing, Xneg)
	m_big = mobile_labor_model(d, _pr_demand(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
		financing = TaxFinanced(_PR_G), pricing = MarkupPricing(20.0))
	Xμ = copy(X); Xμ[5] = 0.5 * d.λ[2]                  # μ_2 = 1 + 20·(0.5 − 1) < 0
	@test_throws PricingFeasibilityError OpenCGE._assert_pricing_feasible(m_big, m_big.pricing, Xμ)
	@test occursin("injections", sprint(showerror, PricingFeasibilityError("injections")))
end

# ── Route 2, Stage B: distributed profit (ADR-0030) ──────────────────────────

@testset "pricing: Stage B nests Stage A (s_π = 1) and the baseline (κ = 0)" begin
	d = split_bill(pricing_fixture())
	for (cl, η, fin) in _PR_CELLS
		mA = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(0.5))
		mB = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(0.5; s_pi = 1.0))
		@test bit_identical(mA, solve(mA), mB, solve(mB))
		m0 = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl, financing = fin)
		mκ = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(0.0; s_pi = 0.5))
		@test bit_identical(m0, solve(m0), mκ, solve(mκ))
	end
	# The keyword: zero profit income is not added, nonzero is added untaxed.
	m = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, 1.0; financing = TaxFinanced(_PR_G))
	p = [1.0, 1.1, 0.9]
	E0 = household_expenditure(m.financing, m, 0.8, p, 0.8; external_transfer = 0.01)
	@test household_expenditure(m.financing, m, 0.8, p, 0.8; external_transfer = 0.01,
		profit_income = 0.0) === E0
	@test household_expenditure(m.financing, m, 0.8, p, 0.8; external_transfer = 0.01,
		profit_income = 0.05) ≈ E0 + 0.05 atol = 1e-15
end

@testset "pricing: Stage B budget, account and diagnostics" begin
	d = split_bill(pricing_fixture())
	for (cl, η, fin) in _PR_CELLS, s_pi in (0.0, 0.5)
		m = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, η; closure = cl,
			financing = fin, pricing = MarkupPricing(0.5; s_pi = s_pi))
		s = solve(m)
		@test maximum(abs, equilibrium_residuals(s)) ≤ 1e-9
		X = _pr_vec(m, s)
		can = external_balance_canary(m, X)
		gc = gdp_components(m, s)
		Π = gc.profit
		@test abs(Π) > 1e-3
		# Budget: E = (1 − τ) w ΣL + F + (1 − s_π) Π, the tax on wages only.
		w, blocks = OpenCGE._solution_blocks(m, s)
		L = sum(blocks.L_i)
		wb = gc.wage_bill
		τ = tau_rate(fin, m, s.prices_raw, wb / L, L)
		@test blocks.E ≈ (1 - τ) * wb + s.external_transfer + (1 - s_pi) * Π atol = 1e-12
		# Reported consumption is the augmented household's.
		@test dot(s.prices_raw, s.consumption) ≈ (1 - d.saving_rate) * blocks.E atol = 1e-12
		# Account: saving on the augmented E plus the retained s_π Π; the
		# distributed rent's import content is inside M. The identity closes.
		@test can.S ≈ d.saving_rate * blocks.E atol = 1e-14
		@test can.retained_profit ≈ s_pi * Π atol = 1e-14
		@test retained_profit(m, s) ≈ s_pi * Π atol = 1e-14
		@test abs(can.diff) ≤ 1e-9
		@test gc.wedge ≈ -can.diff atol = 1e-12
		@test profit_share(m, s) ≈ Π / (wb + Π) atol = 1e-14
	end
	# real_income_index: 1 at the base; the real wage-bill index without
	# distributed profit; labour plus (1 − s_π) Π otherwise.
	mb = mobile_labor_model(d, Shocks(ones(3), ones(3), zeros(3)), 0.5, 0.5, 0.9, 1.0;
		closure = :fixed, financing = TaxFinanced(_PR_G))
	base = solve(mb)
	@test real_income_index(base, base) ≈ 1.0 atol = 1e-15
	@test profit_share(mb, base) == 0.0
	for s_pi in (1.0, 0.5)
		m = mobile_labor_model(d, _pr_supply(), 0.5, 0.5, 0.9, 1.0; closure = :fixed,
			financing = TaxFinanced(_PR_G), pricing = MarkupPricing(0.5; s_pi = s_pi))
		s = solve(m)
		gc, gc0 = gdp_components(m, s), gdp_components(mb, base)
		inc = gc.wage_bill + (1 - s_pi) * gc.profit
		@test real_income_index(s, base) ≈
			(inc / gc0.wage_bill) / (OpenCGE.cpi(s) / OpenCGE.cpi(base)) atol = 1e-14
	end
end
