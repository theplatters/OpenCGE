# STATUS.md — OpenCGE status (2026-10-15)

Plain prose, deliberately short; this is not a generated status board. The
decision record is `DECISIONS.md`; the paper-side evidence (run manifests,
closure catalogue, ADR bodies) lives in the BeyondHulten paper repository.

## What this is

OpenCGE is the production-network CGE kernel of the BeyondHulten project in
the Baqaee–Farhi lineage: a static N-sector economy with CES
production/network structure, pluggable labour-market, financing and pricing
closures, an explicit external account, and national-accounts (Törnqvist)
measurement. It was split out of the BeyondHulten monorepo in October 2026
with its git history preserved; the paper, registry and run records stayed
there.

## Closure landscape

(source of truth: the paper repository's closure catalogue)

- Labour axis: `BF`, `BETA`, `DELTA` are `tested`; `ALPHA` and `GAMMA` are
  `implemented`; `ZETA` (one-sided wage floor with rationed employment) is an
  `idea` and has no constructor (`closure_constructor(:ZETA) === nothing`).
- Financing axis: `F1` (preference reallocation), `F2` (tax financed) and `F3`
  (external debt) are all `tested`.
- Pricing axis (v0.2.0, ADR-0029/ADR-0030): `P1` (`CapacityPricing`,
  capacity / Verdoorn) and `P2` (`MarkupPricing`, retained or distributed
  markup) are in the kernel with synthetic contract tests; they have no
  catalogue entry and no executed generation yet.

## What is validated

The executed evidence is the 5×3 matrix (five labour × three financing cells)
on the full-71 A-bill calibration, run in ten generations; the citable
generation is **`matrix_5x3_v10`** (33/33 cells, all gates pass), plus the
supply-side identification arm `supply_etas_s1`/`_prog`/`_unif` (189 cells,
the η_s recovery design). The kernel behaviours those runs rest on are pinned
by this repository's synthetic test suite: all-N goods-market clearing with
the endogenous external transfer `F`, the η = 0 sectoral-wage endpoint,
financing neutrality of F2/F3, real-wage labour supply, verified scale
determinacy of the fixed-wage system, the monotone ~1e-10 polish,
CES-consistent leak valuation, and the external-account acceptance identity
`S + T_int + M − (I+X) = F + B_gov`.

The pricing closures are pinned on synthetic fixtures only
(`tests/test_pricing_closure.jl`): exact nesting (`δ = 0`, `κ = 0`,
`s_pi = 1` bit-identical), the Cobb–Douglas price identity, the labour FOC at
marginal cost, and the external-account identity with the rent booked, on a
fixture with a hand-split intermediate bill so the bill-factor sites are
exercised. No real-table pricing cell has been run.

## Open items

- The recombination of the two labour margins (allocation η × supply
  elasticity η_s) is a **deliberate omission**, recorded in the paper
  repository's closure catalogue.
- The source IO table's own production-vs-expenditure residual (5.387 % of
  GDP) is disclosed, not modelled.
- The government-side recycling of `T_int` (product taxes on intermediate use)
  is unbooked.
- ADR-0021 (wage structure in the fixed-wage closure) is still `proposed`.
- Pricing door, deferred by ADR-0029: the `μ ≥ 1` domain question (the
  symmetric markup gives contracting sectors `μ < 1`, so the "rent" can be
  negative), the mobile-gauge ladder for both routes, and the real-table
  check of the bill-factor sites (paper repository). `pricing_admissibility`
  reports multi-start spread and Jacobian conditioning; the gate on it lives
  with the paper repository. The ADR-0030 body is not written yet.
