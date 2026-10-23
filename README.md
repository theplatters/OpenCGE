# OpenCGE

A production-network CGE model in the Baqaee–Farhi lineage: N-sector economy,
CES production/network structure, pluggable labour-market, financing and
pricing closures, an explicit external account, and national-accounts
(Törnqvist) measurement of real GDP.

OpenCGE is the kernel of the [BeyondHulten](https://github.com/theplatters/BeyondHulten)
project on green public investment, split out into its own package so that the
model can be developed, tested and cited independently of the paper, its data
and its run records.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/theplatters/OpenCGE", rev = "v0.4.1")
```

For development, clone and `dev` the checkout:

```julia
Pkg.develop(path = "/path/to/opencge")
```

Requires Julia ≥ 1.9. Optional plotting (the GLMakie extension) loads
lazily: `Pkg.add("GLMakie")` and `using GLMakie` in the same session.

## Quick start

```julia
using OpenCGE, DataFrames, LinearAlgebra

# A minimal two-sector calibration (the same fixture the test suite uses).
io = DataFrame("Sektoren" => ["a", "b"],
    "Letzte Verwendung von Gütern zusammen" => [1.0, 1.0])
data = Data(io, Matrix{Float64}(I, 2, 2), [.5, .5], [.5, .5], [1., 1.],
    [.5, .5], [.5, .5], [1., 1.], [0.5, 0.5])

shocks = Shocks(ones(2), ones(2), zeros(2))          # demand, supply, manna
model  = mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0)   # θ, ϵ, σ, η
sol    = solve(model)
sol.prices_raw, sol.quantities, sol.external_transfer

# The closure catalogue: ids, their axis, and the constructor implementing them.
closure_ids()                  # [:BF, :ALPHA, :BETA, :GAMMA, :DELTA, :ZETA, :F1, :F2, :F3, :P1, :P2]
closure_axis(:BETA)            # :labor
closure_constructor(:BETA)     # beta_model

# A pricing closure (ADR-0029, ADR-0031): demand enters the price block;
# the capacity exponent is δ scaled by each sector's value-added share.
solve(mobile_labor_model(data, shocks, 0.5, 0.5, 0.9, 1.0; pricing = CapacityPricing(0.5)))

# First-order welfare measures of a solved cell against a base (read-only).
base = solve(mobile_labor_model(data, Shocks(ones(2), ones(2), zeros(2)), 0.5, 0.5, 0.9, 1.0))
real_consumption_ev(sol, base)                 # money-metric household read, E / π(p)
real_final_demand_ces(sol, base; lambda = 1.0) # CES index of the c + λg bundle
welfare_leisure(sol, base).welfare             # bundle index less w·ΔN / V₀
```

`gamma_model(...)`, `delta_model(...)`, and `closure=:fixed` enforce a fixed
real wage, `w / cpi(sol) = 1`, with endogenous employment and household external
transfer `F = 0`. The solver retains every goods-market and unit-cost equation;
if it cannot also satisfy the real-wage restriction it raises
`FixedRealWageFeasibilityError`. Pricing and productivity shocks can make these
conditions incompatible. Fixing the nominal wage as the gauge does not remove
this economic restriction.

For `:fixed`, `equilibrium_residuals` now appends the real-wage residual after
the N cost and N goods residuals (2N+1 residuals; the unknown vector remains
`[p; y]`). The solver accepts only a common root of all equations.

To reproduce the earlier wage-normalized pricing formulation, explicitly use
`closure=NominalWageClosure()` (or `:nominal`) in `mobile_labor_model`. It fixes
`w = 1` and `F = 0` while allowing the real wage to move. Earlier GAMMA/DELTA
pricing results need new runs under the fixed-real-wage formulation; some cells
may have no compatible equilibrium. `delta_model` accepts a `pricing` keyword.

`read_data("I-O_DE2019_formatiert.csv"; datadir = "...")` loads a real
input–output table (e.g. the German 2019 use table) into the same `Data`
type. The package itself needs **no data**; `read_data` is only exercised at
runtime by whoever supplies the table.

## Repository layout

| Path | What it is |
| --- | --- |
| `src/OpenCGE.jl` | Package entry point (module `OpenCGE`, the exports). |
| `src/core/` | Accounting, calibration, technology, equilibrium, diagnostics. |
| `src/closures/` | Labour, financing and pricing closure plugins, plus the status-free id→constructor map `registry.jl`. |
| `ext/` | `OpenCGEGLMakieExt` (GLMakie plotting, loaded only when GLMakie is loaded) and the matrix-plot helpers. |
| `tests/` | The test suite (synthetic fixtures only) and its entry `tests/runtests.jl`; `test/runtests.jl` is the `Pkg.test` shim. |
| `DECISIONS.md` | The decision record (ADR/DE/IDEA one-liners). |
| `STATUS.md` | Short plain-prose status. |
| `AGENTS.md` | Operating contract for agents and contributors. |

## Comparing closures

`plot_closure_comparison` uses an existing `MatrixDataset` of relative
sectoral results; it does not run models. Select sectors in the desired order
and either labour ids (expanded across F1–F3) or exact run ids:

```julia
using OpenCGE, GLMakie

# ds is a MatrixDataset constructed from your solved cells and baseline.
sectors = [3, 1, 2]                     # indices, sector ids, or display names
closures = [:BF, :BETA]                 # all three financing cells for each
fig = plot_closure_comparison(ds, sectors, closures;
    sector_labels = ["Construction", "Machinery", "Electrical equipment"])
save_closure_comparison("closures.png", ds, sectors, closures;
    sector_labels = ["Construction", "Machinery", "Electrical equipment"])

# Exact run ids select just those cells, in the requested order:
fig = plot_closure_comparison(ds, sectors, [ds.cells[1].run_id, ds.cells[2].run_id])
```

Each sector has output, household consumption (W_D), and price bands, using
the same orange, green, and purple as the existing panel figures. Translucent
bands span the exact minimum and maximum across the selected cells, with
crisp endpoint caps even for unchanged responses. Marker shapes and
filled/outlined styles identify closures in a compact legend. Only nearby
markers (within 1.5% of the displayed x-axis span) receive small symmetric
vertical offsets inside the band; their x-values are preserved.
These ranges compare closures, not confidence intervals or
elasticity variants. Missing, ambiguous, or undefined selected results are
errors. `closure_comparison_data` exposes the changes and extrema headlessly.

Defaults are 17 cm width and 10-point text at that width. The PNG helper
exports at 600 dpi and records the physical resolution; embed it at 17 cm
to retain the intended font size. `width_cm`, `fontsize`, and `dpi` are
configurable (minimum 8-point text and 300 dpi). Bands stay 9 points thick
and sector rows stay 56 points high as closures are added; extra closures
add only compact legend rows. Long sector labels can increase row height.
Split long sector lists across figures to fit an A4 page without
reducing the text. Sector names wrap without truncation. No title or run
metadata is added to the publication figure.

## Comparing closures in an output/price scatter

`plot_closure_scatter` compares explicit labour–financing closure pairs on the
matrix panel's output-change/price-change axes. All sectors are circular points.
Each pair has its own color, and point area is proportional to pre-shock Domar
weight. A sector keeps the same size across every pair. Shocked sectors use
darker opaque shades; other sectors use lighter, more transparent shades.

```julia
using OpenCGE, GLMakie

# ds is an existing MatrixDataset; programme follows its sector order.
pairs = [(:BF, :F1), (:BETA, :F2), (:DELTA, :F3)]
fig = plot_closure_scatter(ds, pairs; programme = programme)
save("closure_scatter.png", fig; px_per_unit = 2)

# Attach names directly to the pairs for the top legend.
named_pairs = ["Benchmark" => (:BF, :F1), "Alternative" => (:BETA, :F2)]
fig = plot_closure_scatter(ds, named_pairs; programme = programme)

# Optionally supply the calibrated pre-shock Domar weights explicitly.
fig = plot_closure_scatter(ds, pairs;
    programme = programme, domar_weights = data.λ,
    unshocked_alpha = 0.25, max_markersize = 24)
```

Only the supplied pairs are selected, in order; `:BF => :F1` pairs are also
accepted. Missing or ambiguous pairs are errors. Default point areas use
baseline `prices .* quantities`, proportional to Domar weights because the
common baseline GDP denominator cancels. Sizes depend on pre-shock values,
including for unshocked sectors. `max_markersize` sets the largest diameter;
all others scale by the square root of relative weight, so areas carry the
correct ratios. Every nonzero programme entry marks a shocked sector.
`closure_colors` optionally provides one distinct base color per pair.
Named and unnamed pairs may be mixed; unnamed pairs keep their closure ids.
`closure_labels` supplies or overrides all top-legend names in selection order.
`closure_scatter_data` exposes percentage changes, the shock mask, and sizing
weights without GLMakie.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite is deliberately synthetic-only: no `data/`, `runs/`, registry or
script machinery is touched, so it is green on a bare checkout. The
real-table regression tests (calibration anchors, external-account cells,
measurement anchors) live in the paper repository, which holds the data and
the executed run manifests.

## Relationship to the paper repository

[github.com/theplatters/BeyondHulten](https://github.com/theplatters/BeyondHulten)
is the paper repository: manuscript, preregistered experiment designs, run
manifests, the closure status catalogue and the real-data test suite. It
depends on OpenCGE for the model; this repository keeps only the kernel, its
tests and its decision record. Full ADR bodies for the decisions summarised in
`DECISIONS.md` are either in this repository's git history or in the paper
repository, as each entry states.

## License

MIT — see `LICENSE`.
