# AGENTS.md — operating contract

## What this repository is

OpenCGE: the production-network CGE kernel (Baqaee–Farhi lineage) of the
BeyondHulten project. One kernel, one package, synthetic-only tests. There is
deliberately **no registry/status-board machinery** here — the paper
repository (https://github.com/theplatters/BeyondHulten) carries the closure
status catalogue, run manifests and the real-data test suite.

## Rules

- **One kernel.** Everything model-related lives in `src/` and is reachable
  from `src/OpenCGE.jl`. No second model implementation, no copied residual
  systems, no zombie files. Closure changes go into `src/closures/` and are
  exposed through the id→constructor map (`closure_ids`, `closure_axis`,
  `closure_constructor`).
- **Tests are synthetic-only** and live in `tests/` (entry `tests/runtests.jl`,
  `Pkg.test` shim in `test/runtests.jl`). A kernel test must never read
  `data/`, `runs/`, a registry, scripts or experiment designs, and must not
  `include()` anything outside `tests/`. Real-table regression tests belong in
  the paper repository.
- **Public API is stable.** The paper repository depends on the exports of
  `OpenCGE`; do not remove or silently change them. Additive changes are fine.
- **Decisions are recorded in `DECISIONS.md`.** Append a new entry (with the
  next id) when a decision changes; never rewrite old entries — supersession
  is stated in the new entry. Dead ends and untested ideas get their line
  there too. No plan files in the repository.
- **Statuses are not restated in code or comments.** If something needs a
  status, it belongs to the paper repository's catalogue.
- **Commit conventions.** `closure(<id>): …` for closure formulation or
  implementation changes, `exp(<design>): run <id>` for executed runs (the
  paper repository), `adr(NNNN): …` / `deadend(NNNN): …` / `docs: …` for
  tracking records; otherwise short imperative subjects. Keep commits
  focused; never mix generated artifacts with model changes.
- **Frozen history is read-only.** The pre-split history documents the
  original BeyondHulten layout; do not rewrite it.
- **House rules.** Empty-vector broadcasting is a trap: after NonlinearSolve
  loads, `Int.(Any[]) === Any[]` — use `Vector{Int}(x)` for values that can be
  empty. `final_demand_split` returns sector × category matrices:
  `vec(sum(fd.tot; dims = 1))` is a category vector. Scale determinacy is a
  property of a fixture: closed fixtures (m = s = 0) are rejected at the
  fixed-wage η = 1 system; give a fixture `s > 0` when a test needs a
  determinate fixed-wage solve. A synthetic fixture cannot test the
  intermediate-bill identity (it is consistent under any split of the bill) —
  that check lives in the paper repository's real-table suite.

## Commands

| Purpose | Command |
| --- | --- |
| Install dependencies | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Load check | `julia --project=. -e 'using OpenCGE'` |
| Test suite | `julia --project=. -e 'using Pkg; Pkg.test()'` |

`Pkg.test()` must be green before any commit lands. Never edit `src/` while
`Pkg.test` is precompiling — the failure is unattributable.
