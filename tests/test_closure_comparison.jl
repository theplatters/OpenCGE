using OpenCGE
using Test

"""Range signs and ordering deliberately differ across the financing cells."""
function closure_comparison_fixture()
    sectors = ["s1", "s2", "s3"]
    labels = ["Construction", "Electrical equipment", "Transport equipment"]
    baseline = MatrixBaseline(ones(3), ones(3), ones(3), ones(3), zeros(3), zeros(3))
    cells = MatrixCellData[]
    for (lab, shift) in (("BF", 0.0), ("BETA", 0.25)), (k, fin) in enumerate(("F1", "F2", "F3"))
        output = [2.0 + k, -5.0 + k, -2.0 + k] .+ shift
        consumption = [-2.0 + k, -3.0 + k, 1.0 + k] .+ shift
        price = [0.0, 0.1 * (k - 2), 0.0]
        push!(cells, MatrixCellData("$lab-$fin", lab, fin, 1.0, 0.5, "ok",
            copy(sectors), copy(labels), 1 .+ price ./ 100, 1 .+ output ./ 100,
            ones(3), 1 .+ consumption ./ 100, ones(3), ones(3),
            Dict{String,Any}(), Dict{String,Any}()))
    end
    return MatrixDataset("synthetic-comparison", baseline, cells,
        Dict{String,Any}(), Dict{String,String}())
end

@testset "Closure-comparison data (headless)" begin
    ds = closure_comparison_fixture()
    data = closure_comparison_data(ds, [3, 1, 2], [:BF])
    @test data.sectors == ["s3", "s1", "s2"]
    @test data.labels == ["Transport equipment", "Construction", "Electrical equipment"]
    @test data.run_ids == ["BF-F1", "BF-F2", "BF-F3"]
    @test data.variables == (:output, :consumption, :price)
    @test size(data.changes) == (3, 3, 3)
    @test data.minima[:, 1] ≈ [-1, 3, -4]
    @test data.maxima[:, 1] ≈ [1, 5, -2]
    @test data.minima[2, 2] ≈ -1
    @test data.maxima[2, 2] ≈ 1
    @test data.minima[2, 3] == data.maxima[2, 3] == 0
    @test closure_comparison_data(ds, ["Electrical equipment", "s1"], ["BF-F3", "BETA-F1"]).run_ids ==
        ["BF-F3", "BETA-F1"]
    data = closure_comparison_data(ds, [1], [:BETA, :BF]; financing = [:F3, :F1])
    @test data.run_ids == ["BETA-F3", "BETA-F1", "BF-F3", "BF-F1"]
    @test data.minima[1, 1] ≈ 3
    @test data.maxima[1, 1] ≈ 5.25
    @test closure_comparison_data(ds, [1], ["BF-F1"]).minima ==
        closure_comparison_data(ds, [1], ["BF-F1"]).maxima
    data = closure_comparison_data(ds, [1], [:BF];
        sector_labels = ["Building"], closure_labels = ["Tax", "Debt", "Preferences"])
    @test data.labels == ["Building"]
    @test data.closure_labels == ["Tax", "Debt", "Preferences"]
    @test_throws DimensionMismatch closure_comparison_data(ds, [1], [:BF]; sector_labels = ["a", "b"])
    @test_throws DimensionMismatch closure_comparison_data(ds, [1], [:BF]; closure_labels = ["a"])
    for (sectors, closures) in ((Int[], [:BF]), ([1], Symbol[]), ([0], [:BF]),
            ([4], [:BF]), (["unknown"], [:BF]), ([1, "s1"], [:BF]),
            ([1], [:unknown]), ([1], [:BF, "BF-F1"]))
        @test_throws ArgumentError closure_comparison_data(ds, sectors, closures)
    end
    @test_throws ArgumentError closure_comparison_data(ds, [1], [:BF]; financing = [:F1, :F1])
    @test_throws ArgumentError closure_comparison_data(ds, [1], [:BF]; financing = Symbol[])
    missing = deepcopy(ds); popfirst!(missing.cells)
    @test_throws ArgumentError closure_comparison_data(missing, [1], [:BF])
    ambiguous = deepcopy(ds); push!(ambiguous.cells, first(ambiguous.cells))
    @test_throws ArgumentError closure_comparison_data(ambiguous, [1], [:BF])
    bad = deepcopy(ds); bad.cells[2].sectors[1:2] = ["s2", "s1"]
    @test_throws ArgumentError closure_comparison_data(bad, [1], [:BF])
    bad = deepcopy(ds); pop!(bad.cells[2].prices)
    @test_throws DimensionMismatch closure_comparison_data(bad, [1], [:BF])
    for value in (NaN, Inf)
        bad = deepcopy(ds); bad.cells[2].consumption[1] = value
        @test_throws ArgumentError closure_comparison_data(bad, [1], [:BF])
        @test size(closure_comparison_data(bad, [2], [:BF]).changes) == (1, 3, 3)
    end
    if Base.get_extension(OpenCGE, :OpenCGEGLMakieExt) === nothing
        @test_throws ErrorException plot_closure_comparison(ds, [1], [:BF])
        @test_throws ErrorException save_closure_comparison("unused.png", ds, [1], [:BF])
    end
end

if get(ENV, "OPENCGE_TEST_GLMAKIE", "false") == "true"
    using GLMakie
    GLMakie.activate!(; visible = false)

    @testset "Publication closure-comparison graphics" begin
        ds = closure_comparison_fixture()
        fig = plot_closure_comparison(ds, [3, 1, 2], [:BF])
        ax = only(filter(b -> b isa Axis, fig.content))
        @test size(fig.scene, 1) == round(Int, 17 / 2.54 * 72)
        @test ax.xticklabelsize[] == ax.yticklabelsize[] == ax.xlabelsize[] == 10
        @test ax.yreversed[]
        limits = ax.finallimits[]
        locations = GLMakie.Makie.get_tickvalues(ax.xticks[],
            Float64(limits.origin[1]), Float64(limits.origin[1] + limits.widths[1]))
        @test length(locations) >= 4
        @test any(<(0), locations)
        @test ax.yticks[][2] == ["Transport equipment", "Construction", "Electrical equipment"]
        bars = filter(p -> p isa GLMakie.Makie.BarPlot, ax.scene.plots)
        @test length(bars) == 3
        @test [p.color[] for p in bars] == GLMakie.Makie.wong_colors()[[2, 3, 4]]
        @test bars[1].fillto[] ≈ [-1, 3, -4]
        @test last.(bars[1][1][]) ≈ [1, 5, -2]
        ticks = filter(p -> p isa GLMakie.Makie.LineSegments, ax.scene.plots)
        @test length(ticks) == 9
        @test all(length(p[1][]) == 6 for p in ticks)
        # Flat price responses still have three separate visible closure ticks.
        @test first(ticks[7][1][])[2] != first(ticks[8][1][])[2]
        legends = filter(b -> b isa Legend, fig.content)
        @test length(legends) == 2
        @test all(l.labelsize[] == 10 for l in legends)
        for legend in legends
            box = legend.layoutobservables.computedbbox[]
            @test all(box.origin .>= 0)
            @test all(box.origin .+ box.widths .<= size(fig.scene))
        end
        more = plot_closure_comparison(ds, [1, 2, 3], [:BF, :BETA])
        @test size(more.scene, 2) > size(fig.scene, 2)
        multi_legend = last(filter(b -> b isa Legend, more.content))
        @test multi_legend.orientation[] == :vertical
        @test multi_legend.nbanks[] == 3
        ext = Base.get_extension(OpenCGE, :OpenCGEGLMakieExt)
        @test ext._cc_wrap("Electrical equipment manufacturing", 16) == "Electrical\nequipment\nmanufacturing"
        @test replace(ext._cc_wrap("Supercalifragilisticexpialidocious", 10), "\n" => "") ==
            "Supercalifragilisticexpialidocious"
        tiny = deepcopy(ds)
        for cell in tiny.cells
            cell.quantities .= 1 + 1e-8
            cell.consumption .= 1
            cell.prices .= 1
        end
        smallfig = plot_closure_comparison(tiny, [1], [:BF])
        smallax = only(filter(b -> b isa Axis, smallfig.content))
        @test smallax.xtickformat[] == "{:.2e}%"
        @test smallax.finallimits[].widths[1] < 1e-4
        @test_throws ArgumentError plot_closure_comparison(ds, [1], [:BF]; width_cm = 5)
        @test_throws ArgumentError plot_closure_comparison(ds, [1], [:BF]; fontsize = 6)
        @test_throws ArgumentError save_closure_comparison("unused.png", ds, [1], [:BF]; dpi = 72)
        @test_throws ArgumentError save_closure_comparison("unused.pdf", ds, [1], [:BF])
        mktempdir() do dir
            path = joinpath(dir, "comparison.png")
            @test save_closure_comparison(path, ds, [3, 1, 2], [:BF]; dpi = 300) == path
            image = GLMakie.FileIO.load(path)
            @test size(image, 2) == round(Int, size(fig.scene, 1) * 300 / 72)
            # Inspect the PNG pHYs chunk: physical print size is encoded.
            bytes = read(path)
            chunk = findfirst(codeunits("pHYs"), bytes)
            @test chunk !== nothing
            offset = last(chunk)
            xppm = foldl((a, b) -> (a << 8) | UInt32(b), bytes[offset+1:offset+4]; init = UInt32(0))
            @test xppm == round(UInt32, 300 / 0.0254)
            @test bytes[offset + 9] == 1
        end
    end
end
