using BeyondHulten
using TOML
using CSV
using DataFrames
using Test

# Contract tests for the Phase 3 experiment entry point (ADR-0006) against
# temporary roots: preregistration gate, cell execution, manifest/index/
# scenario lifecycle, and run-dir immutability. All solves use tiny_fixture()
# with an explicit 2-sector smoke programme (no data/ needed, deterministic,
# fast). Parsed TOML structures are the contract (TOML.print may render some
# tables non-inline, so nothing text-matches the manifest). The scenario-row
# String7 failure below was a real bug: update_scenario_row assigned long
# evidence strings into narrow inferred columns, so every --design run on a
# fresh registry failed outside the cell manifest; run.jl now reads the
# registry with stringtype = String.

# run.jl is include-safe (PROGRAM_FILE guard); skip the re-include when another
# test file already loaded it (both are included into Main by runtests.jl).
isdefined(Main, :run_design) || include(joinpath(@__DIR__, "..", "experiments", "run.jl"))

"""Smoke design TOML: explicit 2-sector programme, single-rung reference."""
function smoke_design_toml()::String
    return """
    schema_version = 1
    design = "smoke"
    description = "2-sector smoke design for tooling tests only"
    smoke = true

    [data]
    vintage = "tiny-fixture"
    drops = []
    calibration_root = "cbase2"
    calibration_artifacts = []

    [programme]
    explicit = [0.7, 0.3]
    total_eur_m = 1.0
    f1_shift = "tilt_g0_over_c0"

    [reference]
    labor = "BF"
    eta = 0.0
    theta = 1.0
    epsilon = 0.5
    sigma = 0.9
    exo_scale_steps = 0
    thetas = [1.0]

    [gates]
    residual_tol = 1e-6
    budget_tol = 1e-9
    labour_tol = 1e-6
    wage_tol = 1e-8

    [cells.smoke-BF-F1]
    labor = "BF"
    financing = "F1"
    eta = 0.0
    eta_s = 0.0
    theta = 1.0
    epsilon = 0.5
    sigma = 0.9
    note = "positive smoke cell"

    [cells.smoke-DELTA-F1]
    labor = "DELTA"
    financing = "F1"
    eta = 1.0
    eta_s = 0.0
    theta = 1e-4
    epsilon = 1e-4
    sigma = 1e-4
    delta_epsilon = 1e-4
    note = "expected scale-indeterminacy failure"
    """
end

"""Fresh temp root with a minimal registry/scenarios.csv and the smoke design."""
function smoke_root()::String
    tmp = mktempdir("/tmp/opencode")
    mkpath(joinpath(tmp, "registry"))
    mkpath(joinpath(tmp, "experiments", "designs"))
    write(joinpath(tmp, "experiments", "designs", "smoke.toml"), smoke_design_toml())
    scen = DataFrame(
        run_id = ["smoke-BF-F1", "smoke-DELTA-F1"],
        design = ["smoke", "smoke"],
        status = ["planned", "planned"],
        labor = ["BF", "DELTA"],
        financing = ["F1", "F1"],
        eta = ["TBD", "TBD"], eta_s = ["TBD", "TBD"],
        theta = ["TBD", "TBD"], epsilon = ["TBD", "TBD"], sigma = ["TBD", "TBD"],
        shock = ["impulses.csv", "impulses.csv"],
        magnitude = ["1.0", "1.0"],
        data_vintage = ["tiny", "tiny"],
        evidence = ["", ""], commit = ["", ""],
        actor = ["test", "test"], notes = ["", ""])
    CSV.write(joinpath(tmp, "registry", "scenarios.csv"), scen)
    return tmp
end

"""Read scenarios.csv rows keyed by run_id (all values as strings)."""
function smoke_scenarios(root::AbstractString)::Dict{String,Dict{String,String}}
    raw = DataFrame(CSV.File(joinpath(root, "registry", "scenarios.csv");
        stringtype = String))
    return Dict(string(r.run_id) =>
        Dict(string(c) => string(coalesce(r[c], "")) for c in names(raw)) for r in eachrow(raw))
end

@testset "run manifests: preregistration gate refuses before creating run dirs" begin
    root = smoke_root()
    runs_dir = joinpath(root, "runs")
    # Unpreregistered design: refuses without touching the runs directory.
    err = try
        run_design("smoke"; root = root, runs_dir = runs_dir, data = tiny_fixture())
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("not preregistered", sprint(showerror, err))
    @test !isdir(runs_dir)
    # Pinned design, then tampered: refuses again, still no run dir.
    preregister_design("smoke"; root = root, actor = "test")
    dpath = joinpath(root, "experiments", "designs", "smoke.toml")
    bytes = read(dpath)
    write(dpath, vcat(bytes, Vector{UInt8}("# tamper\n")))
    err2 = try
        run_design("smoke"; root = root, runs_dir = runs_dir,
            cell = "smoke-BF-F1", data = tiny_fixture())
        nothing
    catch e
        e
    end
    @test err2 isa ErrorException
    @test occursin("changed since preregistration", sprint(showerror, err2))
    @test !ispath(joinpath(runs_dir, "smoke-BF-F1"))
    # Restored byte-identical: the gate accepts the matching record.
    write(dpath, bytes)
    ok, _ = preregistration_status("smoke"; root = root)
    @test ok
end

@testset "run manifests: positive cell writes a schema-valid manifest" begin
    root = smoke_root()
    runs_dir = joinpath(root, "runs")
    preregister_design("smoke"; root = root, actor = "test")
    res = run_design("smoke"; root = root, runs_dir = runs_dir,
        cell = "smoke-BF-F1", data = tiny_fixture(), actor = "test")
    @test res == Dict("smoke-BF-F1" => "executed")
    rundir = joinpath(runs_dir, "smoke-BF-F1")
    man = TOML.parsefile(joinpath(rundir, "manifest.toml"))
    @test man["schema_version"] == 1
    @test man["run_id"] == "smoke-BF-F1"
    @test man["design"] == "smoke"
    @test man["cell"] == "smoke-BF-F1"
    @test man["status"] == "executed"
    @test !isempty(man["date"])
    @test man["actor"] == "test"
    prov = man["provenance"]
    @test prov["seed"] == 1234
    @test prov["design_sha256"] == design_sha256("smoke"; root = root)
    for k in ("git_commit", "git_dirty", "julia_version", "manifest_sha256", "data_sha256")
        @test haskey(prov, k)
    end
    scen = man["scenario"]
    @test scen["labour"] == "BF"
    @test scen["financing"] == "F1"
    @test scen["eta"] ≈ 0.0
    @test scen["shock"] == "impulses.csv"
    @test scen["magnitude"] == 1.0
    for k in ("init", "reference", "algorithm")
        @test haskey(man["solver"], k)
    end
    gates = man["gates"]
    @test gates["overall"] == "pass"
    for (g, tol) in (("residual", 1e-6), ("budget", 1e-9), ("labour", 1e-6))
        @test gates[g]["pass"] == true
        @test gates[g]["tolerance"] ≈ tol
        @test isfinite(gates[g]["value"])
        @test gates[g]["value"] < tol
    end
    for k in ("real_gdp", "real_gdp_ref", "real_gdp_rel", "employment",
            "wage", "nominal_gdp", "max_abs_price_dev")
        @test isfinite(man["metrics"][k])
    end
    for k in ("canary_s", "canary_ixm", "canary_diff", "external_balance", "public_budget")
        @test haskey(man["diagnostics"], k)
    end
    @test man["artifacts"]["log"] == "log.txt"
    @test man["artifacts"]["solution"] == "solution.csv"
    @test isfile(joinpath(rundir, "log.txt"))
    # solution.csv: header plus one row per sector.
    lines = filter(!isempty, split(read(joinpath(rundir, "solution.csv"), String), '\n'))
    @test lines[1] == "sector,price,quantity"
    @test length(lines) == 1 + length(tiny_fixture().factor_share)
    # runs/index.csv: exact header plus one row for the run.
    idxpath = joinpath(runs_dir, "index.csv")
    @test split(read(idxpath, String), '\n')[1] ==
        "run_id,date,design,closures,status,gate_summary,headline_metrics,commit"
    idx = DataFrame(CSV.File(idxpath; stringtype = String))
    @test size(idx, 1) == 1
    @test idx[1, :run_id] == "smoke-BF-F1"
    @test idx[1, :status] == "executed"
    @test idx[1, :design] == "smoke"
    # The scenario row is updated in place: status, pins, evidence, commit.
    row = smoke_scenarios(root)["smoke-BF-F1"]
    @test row["status"] == "executed"
    @test row["labor"] == "BF"
    @test row["financing"] == "F1"
    @test occursin("runs/smoke-BF-F1/manifest.toml", row["evidence"])
    @test !isempty(row["commit"])
end

@testset "run manifests: scale-guard failure is manifest-backed" begin
    root = smoke_root()
    runs_dir = joinpath(root, "runs")
    preregister_design("smoke"; root = root, actor = "test")
    res = run_design("smoke"; root = root, runs_dir = runs_dir,
        cell = "smoke-DELTA-F1", data = tiny_fixture(), actor = "test")
    @test res == Dict("smoke-DELTA-F1" => "failed")
    rundir = joinpath(runs_dir, "smoke-DELTA-F1")
    man = TOML.parsefile(joinpath(rundir, "manifest.toml"))
    @test man["status"] == "failed"
    # Exception-failed manifests carry [gates] overall = "fail" without
    # per-gate values, plus an [error] table; the batch continues.
    @test man["gates"]["overall"] == "fail"
    @test !haskey(man["gates"], "residual")
    @test haskey(man, "error")
    @test occursin("scale-indeterminate", man["error"]["message"])
    @test man["artifacts"]["solution"] == ""
    @test !isfile(joinpath(rundir, "solution.csv"))
    @test isfile(joinpath(rundir, "log.txt"))
    idx = DataFrame(CSV.File(joinpath(runs_dir, "index.csv");
        stringtype = String))
    @test idx[1, :run_id] == "smoke-DELTA-F1"
    @test idx[1, :status] == "failed"
    @test smoke_scenarios(root)["smoke-DELTA-F1"]["status"] == "failed"
end

@testset "run manifests: existing run dirs are not overwritten" begin
    root = smoke_root()
    runs_dir = joinpath(root, "runs")
    preregister_design("smoke"; root = root, actor = "test")
    first = run_design("smoke"; root = root, runs_dir = runs_dir,
        cell = "smoke-BF-F1", data = tiny_fixture(), actor = "test")
    @test first == Dict("smoke-BF-F1" => "executed")
    manpath = joinpath(runs_dir, "smoke-BF-F1", "manifest.toml")
    before = read(manpath)
    second = run_design("smoke"; root = root, runs_dir = runs_dir,
        cell = "smoke-BF-F1", data = tiny_fixture(), actor = "test")
    @test second == Dict("smoke-BF-F1" => "refused")
    @test read(manpath) == before
end
