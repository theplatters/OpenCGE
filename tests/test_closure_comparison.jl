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
    @test data.labours == ["BETA", "BETA", "BF", "BF"]
    @test data.financings == ["F3", "F1", "F3", "F1"]
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
        ds = closure_comparison_fixture(; labours = ["BF", "BETA", "GAMMA"])
        fig = plot_closure_comparison(ds, [3, 1, 2], [:BF, :BETA, :GAMMA])
        axes = filter(b -> b isa Axis, fig.content)
        @test length(axes) == 2
        ax, sector_axis = axes
        @test size(fig.scene, 1) == round(Int, 17 / 2.54 * 72)
        @test ax.xticklabelsize[] == ax.yticklabelsize[] == ax.ylabelsize[] == 10
        @test ax.xticks[][2] == repeat(["BF", "BETA", "GAMMA"], 3)
        @test sector_axis.xticks[][2] == ["Transport equipment", "Construction", "Electrical equipment"]
        @test sector_axis.xaxisposition[] == :top
        limits = ax.finallimits[]
        locations = GLMakie.Makie.get_tickvalues(ax.yticks[],
            Float64(limits.origin[2]), Float64(limits.origin[2] + limits.widths[2]))
        @test length(locations) >= 4
        @test any(<(0), locations)
        # One bar per sector × labour closure × variable: 3 × 3 × 3.
        bars = filter(p -> p isa GLMakie.Makie.BarPlot, ax.scene.plots)
        @test length(bars) == 3
        @test sum(length(p[1][]) for p in bars) == 27
        colors = GLMakie.Makie.wong_colors()[[2, 3, 4]]
        @test [GLMakie.Makie.to_color(p.color[]) for p in bars] ==
            [GLMakie.Makie.to_color((c, 0.28)) for c in colors]
        @test all(p.direction[] == :y for p in bars)
        # Output bars, sector-major then labour: extrema over F1–F3 only.
        @test bars[1].fillto[] ≈ [-1, -0.75, -0.5, 3, 3.25, 3.5, -4, -3.75, -3.5]
        @test last.(bars[1][1][]) ≈ [1, 1.25, 1.5, 5, 5.25, 5.5, -2, -1.75, -1.5]
        xs = first.(bars[1][1][])
        @test issorted(xs)
        @test all(first.(bars[1][1][]) .< first.(bars[2][1][]) .< first.(bars[3][1][]))
        # One line style per financing id; each spans its bar at the cell value.
        lines = filter(p -> p isa GLMakie.Makie.LineSegments, ax.scene.plots)
        @test length(lines) == 3
        linestyles = [GLMakie.Makie.convert_attribute(s, GLMakie.Makie.Key{:linestyle}())
            for s in (:solid, :dash, :dot)]
        @test [p.linestyle[] for p in lines] == linestyles
        @test all(length(p[1][]) == 2 * 27 for p in lines)
        f1 = lines[1][1][]
        @test f1[1] ≈ Point2f(xs[1] - 0.4, -1) && f1[2] ≈ Point2f(xs[1] + 0.4, -1)
        @test all(p[1] == p[2] for line in lines for p in Iterators.partition(last.(line[1][]), 2))
        separators = filter(p -> p isa GLMakie.Makie.VLines, ax.scene.plots)
        @test length(only(separators)[1][]) == 2
        legends = cc_test_legends_fit(fig)
        @test all(l.labelsize[] == 10 for l in legends)
        variable_entries = only(legends[1].entrygroups[])[2]
        @test all(e.attributes[:label][] isa GLMakie.Makie.LaTeXString for e in variable_entries)
        @test [String(e.attributes[:label][]) for e in variable_entries] ==
            [raw"$\Delta y$", raw"$\Delta W_D$", raw"$\Delta p$"]
        @test all(only(e.elements) isa PolyElement for e in variable_entries)
        financing_entries = only(legends[2].entrygroups[])[2]
        @test [e.attributes[:label][] for e in financing_entries] == ["F1", "F2", "F3"]
        @test [only(e.elements).attributes[:linestyle][] for e in financing_entries] ==
            [:solid, :dash, :dot]

        renamed = plot_closure_comparison(ds, [1], [:GAMMA, :BF]; financing = [:F3, :F1],
            labour_labels = ["Gamma", "Baqaee–Farhi"], financing_labels = ["High", "Low"])
        renamed_ax = first(filter(b -> b isa Axis, renamed.content))
        @test renamed_ax.xticks[][2] == ["Gamma", "Baqaee–Farhi"]
        @test [e.attributes[:label][] for e in only(last(cc_test_legends_fit(renamed)).entrygroups[])[2]] ==
            ["High", "Low"]
        @test_throws DimensionMismatch plot_closure_comparison(ds, [1], [:BF]; labour_labels = ["a", "b"])
        @test_throws DimensionMismatch plot_closure_comparison(ds, [1], [:BF]; financing_labels = ["a"])
        # Exact run ids group by labour in order of first appearance.
        mixed = plot_closure_comparison(ds, [1, 2], ["BETA-F2", "BF-F1", "BETA-F3"])
        mixed_ax = first(filter(b -> b isa Axis, mixed.content))
        @test mixed_ax.xticks[][2] == ["BETA", "BF", "BETA", "BF"]
        @test [p.linestyle[] for p in mixed_ax.scene.plots if p isa GLMakie.Makie.LineSegments] ==
            linestyles

        many = closure_comparison_fixture(; labours = ["BF", "BETA", "GAMMA", "DELTA", "ZETA"])
        longfig = plot_closure_comparison(many, [1, 2, 3], [:BF, :BETA, :GAMMA, :DELTA, :ZETA];
            width_cm = 12, sector_labels = ["Capacity constrained sector $k" for k in 1:3])
        cc_test_legends_fit(longfig)
        tall = plot_closure_comparison(ds, [1], [:BF]; axis_height_cm = 10)
        @test size(tall.scene, 2) > size(plot_closure_comparison(ds, [1], [:BF]).scene, 2)
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
        smallax = first(filter(b -> b isa Axis, smallfig.content))
        @test smallax.ytickformat[] == "{:.2e}%"
        @test smallax.finallimits[].widths[2] < 1e-4
        for cell in tiny.cells
            cell.quantities .= 1
        end
        flatfig = plot_closure_comparison(tiny, [1], ["BF-F1"])
        flatax = first(filter(b -> b isa Axis, flatfig.content))
        @test [length(p[1][]) for p in flatax.scene.plots if p isa GLMakie.Makie.LineSegments] == [6]
        cc_test_legends_fit(flatfig)
        @test_throws ArgumentError plot_closure_comparison(ds, [1], [:BF]; width_cm = 5)
        @test_throws ArgumentError plot_closure_comparison(ds, [1], [:BF]; fontsize = 6)
        @test_throws ArgumentError plot_closure_comparison(ds, [1], [:BF]; axis_height_cm = 1)
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
