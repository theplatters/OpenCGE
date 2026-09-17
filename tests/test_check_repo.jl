using BeyondHulten
using TOML
using CSV
using DataFrames
using Test

# Contract tests for scripts/check_repo.jl (ADR-0007), the mandatory
# pre/post-batch gate. Coverage map:
#   - checks 1, 2, 5, 6, 7, 8, 9, 10 pass on the real repo (first testset);
#   - check 3 (closure-subtypes) negative: an unregistered AbstractLaborClosure
#     subtype planted in a temp src/ tree is flagged;
#   - check 4 (kernel-reachability) negative: a zombie src/ file outside the
#     include graph is flagged while a reachable one is clean;
#   - check 5 (runs) negatives: a run dir without manifest.toml, and an
#     index/manifest mismatch, are each flagged;
#   - check 6 (runs-visible) negative: a failed row with no manifest that is
#     NOT pre-manifest history is flagged (an allowlisted cbase2-v3 row is not);
#   - check 9 (preregistration) negative: a stale design hash is flagged;
#   - check 10 (wip-limit) negative: two `running` rows are flagged.
# Checks 1, 2, 7, 8 have no temp negative cases here: they are exercised on
# the real repo (board staleness/warnings, registry bidirectionality, tracked
# artifacts, kernel copies).

# check_repo.jl is include-safe (PROGRAM_FILE guard); skip the re-include when
# already loaded (all test files share Main via runtests.jl). run.jl likewise.
isdefined(Main, :check_repo) || include(joinpath(@__DIR__, "..", "scripts", "check_repo.jl"))
isdefined(Main, :run_design) || include(joinpath(@__DIR__, "..", "experiments", "run.jl"))

"""The real repo root (parent of tests/)."""
check_test_root() = normpath(joinpath(@__DIR__, ".."))

"""Fresh temp dir for gate fixtures (approved scratch space)."""
check_tmp() = mktempdir("/tmp/opencode")

"""Write a minimal registry/scenarios.csv with the given (run_id, status, labor, financing) rows."""
function write_scenarios(root::AbstractString, rows::Vector{Tuple{String,String,String,String}})::Nothing
    mkpath(joinpath(root, "registry"))
    df = DataFrame(run_id = [r[1] for r in rows], design = fill("d", length(rows)),
        status = [r[2] for r in rows], labor = [r[3] for r in rows],
        financing = [r[4] for r in rows], eta = fill("0.5", length(rows)),
        eta_s = fill("0", length(rows)), theta = fill("0.5", length(rows)),
        epsilon = fill("0.5", length(rows)), sigma = fill("0.9", length(rows)),
        shock = fill("impulses.csv", length(rows)), magnitude = fill("1.0", length(rows)),
        data_vintage = fill("tiny", length(rows)), evidence = fill("", length(rows)),
        commit = fill("", length(rows)), actor = fill("test", length(rows)),
        notes = fill("", length(rows)))
    CSV.write(joinpath(root, "registry", "scenarios.csv"), df)
    return nothing
end

@testset "check_repo: real repo is clean" begin
    vs = check_repo(check_test_root())
    @test isempty(vs)
end

@testset "check_repo: unregistered closure subtype is flagged" begin
    tmp = check_tmp()
    mkpath(joinpath(tmp, "registry"))
    cp(joinpath(check_test_root(), "registry", "closures.toml"),
        joinpath(tmp, "registry", "closures.toml"))
    mkpath(joinpath(tmp, "src"))
    write(joinpath(tmp, "src", "fake.jl"),
        "struct FakeClosure <: AbstractLaborClosure end\n")
    vs = check_closure_subtypes(tmp)
    @test length(vs) == 1
    @test occursin("FakeClosure", vs[1].message)
    # A registered symbol is clean.
    write(joinpath(tmp, "src", "fake.jl"),
        "struct FixedWageClosure <: AbstractLaborClosure end\n")
    @test isempty(check_closure_subtypes(tmp))
end

@testset "check_repo: zombie kernel file is flagged" begin
    tmp = check_tmp()
    mkpath(joinpath(tmp, "src"))
    write(joinpath(tmp, "src", "BeyondHulten.jl"), "include(\"a.jl\")\n")
    write(joinpath(tmp, "src", "a.jl"), "# reachable\n")
    write(joinpath(tmp, "src", "zombie.jl"), "# not included\n")
    vs = check_kernel_reachability(tmp)
    @test length(vs) == 1
    @test occursin("zombie.jl", vs[1].message)
    @test !any(occursin("a.jl"), string(v.message) for v in vs)
end

@testset "check_repo: run dir without manifest, and index mismatch" begin
    tmp = check_tmp()
    write_scenarios(tmp, [("r1", "planned", "BF", "F1")])
    mkpath(joinpath(tmp, "runs", "r1"))
    write(joinpath(tmp, "runs", "index.csv"),
        EXPECTED_INDEX_HEADER * "\n")
    vs = check_runs(tmp)
    @test any(occursin("no manifest", v.message) for v in vs)
    # Manifest without an index row is the symmetric violation.
    write(joinpath(tmp, "runs", "r1", "manifest.toml"),
        """
        schema_version = 1
        run_id = "r1"
        design = "d"
        cell = "r1"
        status = "executed"
        date = "2026-09-17"
        actor = "test"
        [provenance]
        [scenario]
        labour = "BF"
        financing = "F1"
        [solver]
        [artifacts]
        """)
    write(joinpath(tmp, "runs", "r1", "log.txt"), "journal\n")
    vs2 = check_runs(tmp)
    @test any(occursin("no runs/index.csv row", v.message) for v in vs2)
end

@testset "check_repo: failed rows stay visible unless pre-manifest history" begin
    tmp = check_tmp()
    write_scenarios(tmp, [
        ("new-fail", "failed", "BF", "F1"),
        ("cbase2-v3-ALPHA-F2-mobile", "provisional", "ALPHA", "F2"),
    ])
    vs = check_failed_visible(tmp)
    @test length(vs) == 1
    @test occursin("new-fail", vs[1].message)
end

@testset "check_repo: stale preregistration hash is flagged" begin
    tmp = check_tmp()
    mkpath(joinpath(tmp, "registry"))
    mkpath(joinpath(tmp, "experiments", "designs"))
    write(joinpath(tmp, "experiments", "designs", "d.toml"), "design = \"d\"\n")
    write(joinpath(tmp, "registry", "preregistration.toml"),
        "schema_version = 1\n[designs.d]\ndesign_sha256 = \"deadbeef\"\n")
    vs = check_preregistration(tmp)
    @test length(vs) == 1
    @test occursin("changed since preregistration", vs[1].message)
    # The matching record is clean.
    preregister_design("d"; root = tmp, actor = "test")
    @test isempty(check_preregistration(tmp))
end

@testset "check_repo: WIP limit" begin
    tmp = check_tmp()
    write_scenarios(tmp, [("a", "running", "BF", "F1"), ("b", "planned", "BF", "F2")])
    @test isempty(check_wip_limit(tmp))
    write_scenarios(tmp, [("a", "running", "BF", "F1"), ("b", "running", "BF", "F2")])
    vs = check_wip_limit(tmp)
    @test length(vs) == 1
    @test occursin("WIP limit", vs[1].message)
end
