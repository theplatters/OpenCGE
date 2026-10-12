# OpenCGE

A production-network CGE model in the Baqaee–Farhi lineage: N-sector economy,
CES production/network structure, pluggable labour-market and financing
closures, an explicit external account, and national-accounts (Törnqvist)
measurement of real GDP.

OpenCGE is the kernel of the [BeyondHulten](https://github.com/theplatters/BeyondHulten)
project on green public investment, split out into its own package so that the
model can be developed, tested and cited independently of the paper, its data
and its run records.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/theplatters/opencge", rev = "v0.1.0")
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
closure_ids()                  # [:BF, :ALPHA, :BETA, :GAMMA, :DELTA, :ZETA, :F1, :F2, :F3]
closure_axis(:BETA)            # :labor
closure_constructor(:BETA)     # beta_model
```

`read_data("I-O_DE2019_formatiert.csv"; datadir = "...")` loads a real
input–output table (e.g. the German 2019 use table) into the same `Data`
type. The package itself needs **no data**; `read_data` is only exercised at
runtime by whoever supplies the table.

## Repository layout

| Path | What it is |
| --- | --- |
| `src/OpenCGE.jl` | Package entry point (module `OpenCGE`, the exports). |
| `src/core/` | Accounting, calibration, technology, equilibrium, diagnostics. |
| `src/closures/` | Labour and financing closure plugins, plus the status-free id→constructor map `registry.jl`. |
| `ext/` | `OpenCGEGLMakieExt` (GLMakie plotting, loaded only when GLMakie is loaded) and the matrix-plot helpers. |
| `tests/` | The test suite (synthetic fixtures only) and its entry `tests/runtests.jl`; `test/runtests.jl` is the `Pkg.test` shim. |
| `DECISIONS.md` | The decision record (ADR/DE/IDEA one-liners). |
| `STATUS.md` | Short plain-prose status. |
| `AGENTS.md` | Operating contract for agents and contributors. |

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
