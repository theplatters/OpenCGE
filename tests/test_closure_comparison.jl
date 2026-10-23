using OpenCGE
using Test

"""Range signs and ordering deliberately differ across the financing cells."""
function closure_comparison_fixture(; labours = ["BF", "BETA"])
    sectors = ["s1", "s2", "s3"]
    labels = ["Construction", "Electrical equipment", "Transport equipment"]
    baseline = MatrixBaseline(ones(3), ones(3), ones(3), ones(3), zeros(3), zeros(3))
    cells = MatrixCellData[]
    for (l, lab) in enumerate(labours), (k, fin) in enumerate(("F1", "F2", "F3"))
        shift = 0.25 * (l - 1)
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

    function cc_test_legends_fit(fig)
        legends = filter(b -> b isa Legend, fig.content)
        @test length(legends) == 2
        for legend in legends
            box = legend.layoutobservables.computedbbox[]
            @test all(box.origin .>= 0)
            @test all(box.origin .+ box.widths .<= size(fig.scene) .+ 1)
        end
        return legends
    end

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
        colors = GLMakie.Makie.wong_colors()[[2, 3, 4]]
        @test [GLMakie.Makie.to_color(p.color[]) for p in bars] ==
            [GLMakie.Makie.to_color((c, 0.28)) for c in colors]
        @test all(p.strokewidth[] == 1.0 for p in bars)
        @test [GLMakie.Makie.to_color(p.strokecolor[]) for p in bars] ==
            [GLMakie.Makie.to_color((c, 0.85)) for c in colors]
        @test bars[1].fillto[] ≈ [-1, 3, -4]
        @test last.(bars[1][1][]) ≈ [1, 5, -2]
        caps = filter(p -> p isa GLMakie.Makie.LineSegments, ax.scene.plots)
        @test length(caps) == 3
        # Two flat price ranges have one cap each; the middle range has two.
        @test [length(p[1][]) for p in caps] == [12, 12, 8]
        @test all(p[1][1] == p[2][1] && p[1][2] != p[2][2]
            for cap in caps for p in Iterators.partition(cap[1][], 2))
        markers = filter(p -> p isa GLMakie.Makie.Scatter, ax.scene.plots)
        @test length(markers) == 9
        @test all(all(==(7.5f0), p.markersize[]) && p.strokewidth[] == 0.9f0 for p in markers)
        # Untied output responses stay on the centreline; flat prices separate.
        @test last.(markers[1][1][]) == last.(markers[2][1][]) == last.(markers[3][1][])
        @test length(unique(first(p[1][])[2] for p in markers[7:9])) == 3
        @test sum(first(p[1][])[2] for p in markers[7:9]) / 3 ≈ 1 + 14 / 56
        separators = filter(p -> p isa GLMakie.Makie.HLines, ax.scene.plots)
        @test only(separators)[1][] == [1.5, 2.5]
        legends = cc_test_legends_fit(fig)
        @test all(l.labelsize[] == 10 for l in legends)
        variable_entries = only(legends[1].entrygroups[])[2]
        @test [e.attributes[:label][] for e in variable_entries] ==
            ["Output", "Consumption (W_D)", "Price"]
        @test all(only(e.elements) isa PolyElement for e in variable_entries)
        closure_entries = only(legends[2].entrygroups[])[2]
        for (entry, marker) in zip(closure_entries, markers[1:3])
            element = only(entry.elements)
            @test element isa MarkerElement
            @test GLMakie.Makie.to_spritemarker(element.attributes[:marker][]) == marker.marker[]
            @test GLMakie.Makie.to_color(element.attributes[:markercolor][]) == marker.color[]
            @test Float32(element.attributes[:markerstrokewidth][]) == marker.strokewidth[]
        end
        more = plot_closure_comparison(ds, [1, 2, 3], [:BF, :BETA])
        @test size(more.scene, 2) > size(fig.scene, 2)
        multi_legend = last(filter(b -> b isa Legend, more.content))
        @test multi_legend.orientation[] == :vertical
        @test multi_legend.nbanks[] == 4
        @test multi_legend.patchsize[] == (14, 12)
        @test multi_legend.rowgap[] == 2
        @test multi_legend.colgap[] == 12
        cc_test_legends_fit(more)

        # Near-coincident responses separate without shifting a third, distant
        # closure or the same closures' untied responses in another sector.
        near = deepcopy(ds)
        near.cells[2].quantities[1] = near.cells[1].quantities[1] + 0.0005
        nearfig = plot_closure_comparison(near, [1, 2, 3], [:BF])
        nearax = only(filter(b -> b isa Axis, nearfig.content))
        nearmarkers = filter(p -> p isa GLMakie.Makie.Scatter, nearax.scene.plots)
        @test first(nearmarkers[1][1][])[2] < first(nearmarkers[3][1][])[2] < first(nearmarkers[2][1][])[2]
        @test [p[1][][2][2] for p in nearmarkers[1:3]] ≈ fill(2 - 14 / 56, 3)
        @test abs(first(nearmarkers[2][1][])[2] - first(nearmarkers[1][1][])[2]) * 56 ≈ 2.2 atol = 1e-5

        many = closure_comparison_fixture(; labours = ["BF", "BETA", "GAMMA", "DELTA", "ZETA"])
        for nc in (9, 15)
            selection = [c.run_id for c in many.cells[1:nc]]
            densefig = plot_closure_comparison(many, [1, 2, 3], selection)
            denseax = only(filter(b -> b isa Axis, densefig.content))
            densebars = filter(p -> p isa GLMakie.Makie.BarPlot, denseax.scene.plots)
            densemarkers = filter(p -> p isa GLMakie.Makie.Scatter, denseax.scene.plots)
            @test size(denseax.scene, 2) ≈ size(ax.scene, 2) ≈ 3 * 56
            @test all(p.width[] * size(denseax.scene, 2) / 3 ≈ 9 for p in densebars)
            @test size(densefig.scene, 2) - size(fig.scene, 2) <= 50
            @test size(densefig.scene, 2) < 320
            # Fractional tick locations must retain their decimal places.
            ticklimits = denseax.finallimits[]
            tickvalues = GLMakie.Makie.get_tickvalues(denseax.xticks[],
                Float64(ticklimits.origin[1]), Float64(ticklimits.origin[1] + ticklimits.widths[1]))
            @test denseax.xtickformat[] == (all(isinteger, tickvalues) ? "{:.0f}%" : "{:.1f}%")
            # All closures share flat prices in sector 1: distinct symmetric
            # positions within the same nine-point band, including 15-way ties.
            price_markers = densemarkers[(2nc + 1):3nc]
            ys = [first(p[1][])[2] for p in price_markers]
            @test length(unique(ys)) == nc
            @test sum(ys) / nc ≈ 1 + 14 / 56
            @test (maximum(ys) - minimum(ys)) * 56 <= 9
            @test length(unique((p.marker[], p.color[]) for p in densemarkers[1:nc])) == nc
            cc_test_legends_fit(densefig)
            longfig = plot_closure_comparison(many, [1, 2, 3], selection;
                width_cm = 12, closure_labels = ["Capacity constrained sectoral closure $k" for k in 1:nc])
            cc_test_legends_fit(longfig)
        end
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
        for cell in tiny.cells
            cell.quantities .= 1
        end
        flatfig = plot_closure_comparison(tiny, [1], ["BF-F1"])
        flatax = only(filter(b -> b isa Axis, flatfig.content))
        @test all(length(p[1][]) == 2 for p in flatax.scene.plots if p isa GLMakie.Makie.LineSegments)
        @test all(length(p[1][]) == 1 for p in flatax.scene.plots if p isa GLMakie.Makie.Scatter)
        cc_test_legends_fit(flatfig)
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
