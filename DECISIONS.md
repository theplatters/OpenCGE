# DECISIONS.md — OpenCGE decision record

This file is the kernel's decision record. OpenCGE deliberately carries **no
registry/status-board machinery**: a decision is recorded here in one or two
sentences and appended to when it changes (old entries are never rewritten).
Ids are kept from the originating BeyondHulten project so that the paper, its
tables and the run manifests can be read against this list. Entries marked
"paper repo" keep their full record there; for the kernel-only decisions the
full record lives on as this repository's git history.

## Decisions

Kernel-relevant decisions whose full record belongs to this repository:

- **ADR-0001 — One kernel, no living copies.** The Julia package is the single
  canonical model kernel; the legacy `cbase2`, `bf_replication` and
  `bf_replication2` trees are frozen read-only snapshots and no residual system
  may be built by copying the kernel. (full record: ADR-0001, retained here as
  history)
- **ADR-0005 — Phase-2 layout and backport.** The kernel is laid out as
  `src/core/{accounting,technology,equilibrium,diagnostics}.jl` with closure
  plugins under `src/closures/` and a status-free id→constructor map in
  `src/closures/registry.jl`; the legacy-kernel diffs are ported onto this
  layout with compatibility guarantees (v3 absorption field names/order,
  legacy error surfaces, Törnqvist default). (full record: ADR-0005, retained
  here as history)
- **ADR-0010 — Port the review fixes; keep the BF endpoints.** Ports the
  retained-sector rebuild (retained_io_table/retained_dataset) and the
  N−1+CPI mobile formulation, keeps only the BF endpoints η ∈ {0, 1}, and
  retires the interpolated 0 < η < 1 cases and the ad hoc allocative-efficiency
  wedge; its N−1 clearing decision is later corrected by ADR-0019. (full
  record: ADR-0010, retained here as history)
- **ADR-0011 — Conform the frozen wrangling notebook.** The frozen
  data-wrangling record is edited (under a re-pinned freeze record) to locate
  the final-demand and import/tax blocks **by label** instead of hard-coded
  table positions and to call the canonical implementation, so the submission
  record cannot silently drift from the kernel. (full record: ADR-0011,
  retained here as history)
- **ADR-0014 — Real-wage supply and a verified scale-determinacy guard.**
  Elastic labour supply is defined on the real wage,
  `L^s = L̄·[(w/P)/(w₀/P₀)]^{η_s}` (DE-0004), and the fixed-wage η = 1 system is
  admitted only under a verified determinacy criterion instead of the old
  additive-anchor heuristic. (full record: ADR-0014, retained here as history)
- **ADR-0015 — Monotone polish to ~1e-10.** The residual-gated
  Levenberg–Marquardt polish triggers below 1e-10, runs at 1e-12 tolerances and
  is monotone (a step is kept only if it strictly improves the residual), so
  cell metrics never depend on the solver's stopping point; the acceptance
  gates stay at 1e-6 (fixed) / 1e-5 (mobile). (full record: ADR-0015, retained
  here as history)
- **ADR-0016 — CES-consistent leak valuation.** The intermediate-bill leaks
  `M_int` and `T_int` are valued with the same CES expenditure factor
  `k = p^ε·a^(ε−1)·P^(1−ε)` as the domestic bill (not frozen at 1), so the
  external-account identity closes outside the baseline gauge — required for
  supply-shock scenarios. (full record: ADR-0016, retained here as history)
- **ADR-0017 — Actual-matrix scale determinacy.** The fixed-wage η = 1
  admissibility guard now tests the actual clearing matrix (`I − G` full rank)
  rather than the closed-form round-gain column sums, which both rejected
  determinate systems and mislocated the import-margin feedback (it enters at
  the spending sector). (full record: ADR-0017, retained here as history)
- **ADR-0019 — All-N clearing with an explicit external account.** Every
  regime enforces **all N** goods-market clearings, and the mobile η = 1 system
  carries the endogenous net external transfer `F` as an unknown (canonical
  vector `[p; y; w; F]`, `E = (1−τ)·w·ΣL + F`), with the acceptance identity
  `S + T_int + M − (I+X) = F + B_gov`; this corrects ADR-0010's N−1 decision
  only. (full record: ADR-0019, retained here as history)
- **ADR-0020 — η = 0 endpoint carries sectoral wages (option C).** At η = 0
  the frozen allocation is cost-minimizing at sector-specific wages `w_i` (a
  3N+1 system solving the per-sector first-order condition), so the identity
  closes and the external position is identified; the retired `F = 0` pin left
  it unidentified. (full record: ADR-0020, retained here as history)
- **ADR-0021 — Wage structure in the fixed-wage closure (proposed).** Would
  make the wage vector an explicit dimension of GAMMA by reusing the ADR-0020
  vector-wage apparatus (today the design sets it to its degenerate all-equal
  value, hard-coded as `w = 1`); status `proposed`, nothing implemented.
  (full record: ADR-0021, retained here as history)
- **ADR-0022 — Sectoral labour markets.** BETA is generalised to N sectoral
  labour markets with a per-sector real-wage supply elasticity (3N+1 system,
  `problem_sectoral`), the "labour door" to demand-sensitive prices under
  demand-only shocks, where every single-wage closure has `p` pinned at one;
  the η_s,i = 0 corner nests the ADR-0020 endpoint exactly. (full record:
  ADR-0022, retained here as history)

Kernel-relevant decisions whose full record stays in the paper repository
(BeyondHulten):

- **ADR-0002 — Closure taxonomy and stable IDs.** BF/ALPHA/BETA/GAMMA/DELTA/
  ZETA label the labour axis, F1/F2/F3 the financing axis; `η` is a geometric
  reallocation friction (never a supply elasticity — that is `η_s`) and the
  labels are project ids, not claims of literature equivalence. (full record:
  ADR-0002, paper repo)
- **ADR-0012 — A-bill calibration.** Intermediate demand is charged with the
  **domestic** bill `A_bill` (raw-table row 73), with imported intermediates
  `M_int` moved out of the demand system into the external account, replacing
  the purchaser-price bill and its household-residual clamp. (full record:
  ADR-0012, paper repo)
- **ADR-0013 — Intermediate-bill tax term.** Product taxes on intermediate use
  `T_int` (row 75) are booked as the third component of the purchaser-price
  bill and an explicit external leak; unbooked, the canary is short by exactly
  `ΣT_int` = 2.5745e-2 of GDP. (full record: ADR-0013, paper repo)
- **ADR-0018 — Real-GDP measurement.** Real GDP is the Törnqvist/Divisia
  quantity index of final demand with the signed national-accounts components
  `V = [C_gross, G+programme, I, X, −M_final, −M_int, −T_int]` and the
  wedge/canary duality `gdp_wedge = −canary.diff`; a welfare index is never
  called GDP without a proof of equivalence. (full record: ADR-0018, paper
  repo)
- **ADR-0024 — Programme total.** The investment programme is pinned at
  G0 = 40,300 EUR m (2019 prices, the deflated 2024–2050 horizon mean), the
  magnitude of every executed matrix cell. (full record: ADR-0024, paper repo)
- **OpenCGE@ADR-0029 — The pricing door (Foundation III).** A third closure
  axis makes unit cost output-dependent: `P1` `CapacityPricing(δ)`,
  `A_eff = A·(y/λ)^(−δ)`, rent-free, and `P2` `MarkupPricing(κ)`, `p = μ·c`
  with `μ = 1 + κ(y/λ − 1)` and the rent retained beside saving (Stage A);
  every FOC site reads the marginal cost `p/μ` and `A_eff` through one hook
  set, `δ = 0` and `κ = 0` reproduce the baseline bit for bit, and the
  ADR-0017 guard runs on a NoPricing clone plus a local Jacobian rank check.
  Proposed in the paper repository; implemented in OpenCGE v0.2.0. (full
  record: ADR-0029, paper repo)
- **OpenCGE@ADR-0030 — Distributed profit income (Stage B).**
  `MarkupPricing(κ; s_pi)` with `s_pi < 1` pays `(1 − s_pi)·Π` to the
  household, `E = (1 − τ)·w·ΣL + F + (1 − s_pi)·Π`, untaxed (v1), and the
  canary books saving `s·E + s_pi·Π`; `s_pi = 1` is Stage A exactly.
  Implemented in OpenCGE v0.2.0 ahead of the ADR body, which the paper
  repository has not written yet. (full record: ADR-0030, paper repo)
- **OpenCGE@ADR-0031 — Value-added-scaled capacity exponent.**
  `CapacityPricing(δ)` now reads `A_eff_i = A_i·(y_i/λ_i)^(−δ_i·v_i)` with
  `v_i = factor_share_i` the value-added share of gross output, so `δ` is the
  elasticity of value-added productivity with respect to output (Verdoorn
  coefficient `β = −δ`) and `v_i` turns it into the gross-output exponent;
  the unscaled ADR-0029 form `(y/λ)^(−δ)` mis-states that exponent and is
  superseded. `δ` may be a per-sector vector (a constant vector is the scalar
  bit for bit); `δ = 0` still nests the baseline exactly. Breaking for P1
  results: the same `δ` now gives a weaker exponent. `pricing_technology`
  takes `v` as a fifth argument. (kernel-only; full record: this
  repository's git history)

- **ADR-0034 — Opt-in print layout for matrix panels.** `plot_matrix_panel`
  and `save_matrix_panels` retain the exploratory layout by default and add
  `style = :manuscript`: a 1100 × 560 logical-pixel layout for about 14 cm
  print width, with larger fonts, legends below the axes explaining the
  baseline and central-cell/σ-variant range, and no title, footnotes or run
  id. This changes presentation only; saved PNGs use two pixels per logical
  unit, and `title = ""` also suppresses the exploratory title.

- **ADR-0035 — Compact numbered manuscript panels (v0.3.2).** Supersedes
  ADR-0034's print presentation: 1100 × 520 logical pixels, horizontal
  programme-sector bars numbered in programme order, matching numbers in
  the scatter or above quantity bars, one shared two-row legend, baseline
  axis labels and integer-percent ticks (one decimal below a 3-point span).
  Callers supply curated `sector_labels`; fallback names are cleaned and
  shortened to 25 characters. Zero response series are omitted; error bars
  retain the central-cell/converged-σ-variant min/max meaning. At most two
  non-programme scatter outliers receive nearby labels, accepted only when
  measured data-space text boxes fit without overlapping markers or labels;
  no leader lines are used. Exploratory output is preserved. Presentation
  only: no solver, closure, calibration or data-layer changes. The original
  print-layout entry is renumbered from ADR-0032 to ADR-0034 to avoid the
  paper repository's ADR-0032 (pricing doors on every labour row) and
  ADR-0033 (financing balancing variables); its historical content is retained.
  Project and citation versions are aligned at 0.3.2; tagging awaits approval.

- **ADR-0036 — Genuine fixed real wage with zero household transfer.**
  GAMMA/DELTA (`FixedWageClosure`, `:fixed`) enforce `w/Π(p) = 1` while
  retaining `F = 0`, all N unit-cost equations, and all N goods-market
  clearings. With the nominal gauge `w = 1`, the extra real-wage restriction
  gives 2N+1 equations in 2N unknowns; least squares is used to find a common
  root, and a nonzero residual minimum raises
  `FixedRealWageFeasibilityError` rather than returning a moving real wage or
  introducing a household transfer. This supersedes the nominal-wage reading
  of GAMMA/DELTA in ADR-0021 and the pricing applications of ADR-0029 through
  ADR-0031: those results require new runs and may be incompatible with the
  fixed real wage. The previous equations remain available explicitly as
  `NominalWageClosure()` / `:nominal` for reproduction, with no new closure id.
  Demand-only NoPricing results continue to satisfy the real-wage restriction.
  A shared CES CPI evaluates the geometric-mean limit at `σ = 1` so the
  restriction cannot disappear through the degenerate `1^Inf` expression.
  Synthetic tests cover employment adjustment, infeasible pricing and
  productivity shocks, and compatible pricing with nonuniform prices and
  zero external-account gaps. A numerical failure to find a common root is
  reported as such, not as a global nonexistence proof.

- **ADR-0037 — Closure-comparison range bars.** Add a headless selector and
  GLMakie recipe over existing `MatrixDataset` cells. Each selected sector
  has output, household consumption (W_D), and price bars in the existing
  panel palette (Wong colours 2–4); each bar spans the exact selected-cell
  minimum and maximum, with individual closure ticks in ordered separate
  lanes to preserve coincident responses. Labour ids expand over explicitly
  selected financing ids; exact run ids select individual cells. Missing,
  ambiguous and nonfinite selected responses fail explicitly. Print defaults
  are 17 cm, 10-point text, and 600-dpi PNGs with resolution metadata; height
  grows with the selection and labels wrap. Presentation only; no solver or
  closure formulation changes.

- **ADR-0038 — Compact closure ranges and markers (v0.4.1).** Supersedes
  ADR-0037's separate closure lanes and tick legend. Translucent 9-point
  min–max bands with bordered endpoint caps sit in 56-point sector rows,
  separated by subtle rules. Closure markers use nine shapes, then repeat
  with white fills; a standard marker legend uses up to four columns
  (three for long labels). Only responses connected by gaps of at most
  1.5% of the displayed x-axis span receive symmetric vertical offsets;
  the step is at most 2.2 points and shrinks for large tie groups to stay
  inside the band. X-values and exact extrema remain unchanged. Measured
  legend layout determines figure height while sector and band heights
  remain independent of closure count; long labels wrap. Synthetic graphics
  checks cover near/exact ties, flat responses, 9–15 closures, marker styles,
  compact physical geometry, legend bounds and PNG resolution. Tick precision
  follows the actual tick locations so fractional steps retain their decimals.
  The headless data contract and public plotting signatures remain unchanged. Release as
  v0.4.1, retaining the existing v0.4.0 tag.

- **ADR-0039 — Closure-comparison output/price scatter.** Add a headless
  `closure_scatter_data` selector and GLMakie `plot_closure_scatter` over
  existing `MatrixDataset` cells, sharing the range plot's closure selection
  rules. Plot all sectors' exact percentage changes in output and price on
  shared axes, with sector colors independent of the closure selection and
  distinct marker shapes for up to 15 closures by default. Every nonzero
  programme entry identifies a shocked sector, including negative shocks.
  Shocked sectors use darker opaque markers; other sectors use lighter,
  smaller translucent markers and are drawn first across all closures.
  Sector legends default to shocked sectors, with all/none options. Flat-price
  responses remain visible without jitter; invalid plotted responses fail
  explicitly, while unplotted consumption is irrelevant. Presentation only;
  no solver or closure formulation changes.

- **ADR-0040 — Domar-sized points and explicit closure pairs.** Supersedes
  ADR-0039's sector colors, closure shapes, and shock-dependent marker sizes.
  The scatter accepts explicit labour–financing pairs, each matching one cell,
  and assigns a distinct color to each pair. Every sector is a circular point
  with area proportional to its pre-shock Domar weight, shared across pairs
  and shock groups. Default weights use baseline gross-output values, whose
  ratios equal Domar-weight ratios; explicit `domar_weights` may be supplied.
  Marker diameters scale with square roots of weights. Shocked sectors use
  darker opaque shades; other sectors use lighter translucent shades drawn
  first. Pair colors remain stable under reordering and subsetting within a
  dataset. Legends show closure pairs and explain shock shading and area.
  This replaces the initial scatter input as explicitly requested; the range
  comparison API and all model formulations are unchanged.

- **ADR-0042 — Names attached to scatter closure pairs.** Add optional
  `"Name" => (labour, financing)` entries to the scatter's pair selection.
  Their names appear in the top legend, with closure ids as the fallback
  for unnamed entries in mixed selections. The existing `closure_labels`
  keyword overrides attached names. Label validation retains nonempty,
  unique display names; pair selection, colors, weights and coordinates
  depend only on the underlying closure pair.

- **ADR-0043 — Restore the fixed nominal wage in GAMMA/DELTA.** Supersedes
  ADR-0036: `FixedWageClosure`, `:fixed`, `gamma_model` and `delta_model` fix
  `w = 1` and `F = 0`, with endogenous employment and a moving CPI and real
  wage. Remove the additional `w/Π(p) = 1` restriction and its overdetermined
  least-squares solve: it prevents price-changing shocks from clearing all
  markets. Restore the square system of N unit-cost and N goods-market
  equations in `[p; y]` and the ordinary nonlinear root solver with residual
  checks and bounded polishing. `equilibrium_residuals` again returns 2N
  entries. Keep `NominalWageClosure()` / `:nominal` as compatibility spellings
  for the same equations, retain the exported `FixedRealWageFeasibilityError`
  type without raising it, and keep the shared CPI's geometric limit at
  `σ = 1`. Synthetic regressions cover price-changing pricing and productivity
  shocks, raw versus CPI-normalized wages, and zero household transfer.

## Dead ends and ideas

- **DE-0001 — Unnormalized household demand shifter.** Multiplying CES weights
  without renormalizing breaks the demand identity (measured `Σ ≈ 1.255`);
  shifters must stay normalized.
- **DE-0002 — Dropping a zero-profit equation.** Closes the equation count
  only at the cost of positive profits (sector-1 residual ≈ 0.39).
- **DE-0003 — η = 10⁶ as a fixed-real-wage approximation.** The high-η limit
  collapses employment (L → 0); it is not a fixed-wage economy.
- **DE-0004 — Nominal wage in labour supply.** Violates numeraire homogeneity;
  supply must be on the real wage `w/P` (the basis of ADR-0014).
- **DE-0005 — Fixed-base quantity sum called real GDP.** Not a quantity index;
  use the shared Törnqvist index.
- **DE-0006 — Renormalized variance shares.** The "88.4 %"/"85.9 %" claims are
  shares of reported main effects, not of total variance; incomplete designs
  must fail loudly.
- **DE-0007 — Baqaee–Farhi (2022) wholesale as the base model.** Rejected as a
  base replacement; its sticky sectoral labour complementarity remains the
  seed of the ZETA idea.
- **DE-0008 — Open-economy calibration without a saving rate.** No interior
  equilibrium (import margin without a saving/export block); the saving rate is
  required.
- **DE-0009 — Orphan `src/` files.** Six plausible but unwired files were
  archived instead of left as zombies: every file must be reachable from the
  module.
- **DE-0010 — Multi-stage solver ladder.** Superseded by the single `solve()`
  path with its residual-gated monotone polish (ADR-0015).
- **DE-0011 — Re-gauging the numeraire to give η_s a channel.** Rejected: the
  system is homogeneous of degree 1 in `(p, w, F)`, so re-gauging moves only
  the CPI pin and creates no identification.
- **DE-0012 — Capacity door as an outer fixed point.** Iterating
  `A ← A·(y/λ)^(−δ)` around the solver (probe17/probe20) diverges at damping
  0.6 for δ = −0.5 under Cobb–Douglas; `A_eff` is solved inside the residual
  instead (ADR-0029 #8).
- **IDEA-0001 — Productivity shocks framed as climate change.** An
  interpretation/identification route for the supply-shock arm; unexplored.
- **IDEA-0002 — Five doors to demand-sensitive prices (Verdoorn arm).** A
  ranked menu of mechanisms; only the labour door is implemented (ADR-0022).
- **IDEA-0003 — Sobol decomposition of the sectoral labour response.**
  Variance decomposition of sectoral outcomes over the closure parameters;
  scoped, nothing measured.
- **IDEA-0004 — Sectoral effects without a Sobol.** Cheaper
  affected-vs-unaffected incidence measures and revived sectoral figures;
  scoped, nothing measured.
