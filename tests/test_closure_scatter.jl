using OpenCGE
using Test

function closure_scatter_fixture(; labours = ["BF", "BETA"])
    sectors = ["s1", "s2", "s3", "s4"]
    labels = ["Construction", "Equipment", "Energy", "Services"]
    baseline = MatrixBaseline([1.0, 2.0, 1.0, 1.0], [0.04, 0.16, 0.01, 0.09],
        ones(4), ones(4), zeros(4), zeros(4))
    cells = MatrixCellData[]
    for (l, lab) in enumerate(labours), (f, fin) in enumerate(("F1", "F2", "F3"))
        dq = [2f, -f, 0.2f, -0.1f] .+ (l - 1) / 4
        dp = lab == "BF" ? zeros(4) : [0.2f, -0.1f, 0.01f, -0.02f]
        push!(cells, MatrixCellData("$lab-$fin", lab, fin, 1.0, 0.5, "ok",
            copy(sectors), copy(labels), 1 .+ dp ./ 100, 1 .+ dq ./ 100,
            ones(4), ones(4), ones(4), ones(4), Dict{String,Any}(), Dict{String,Any}()))
    end
    return MatrixDataset("synthetic-scatter", baseline, cells,
        Dict{String,Any}(), Dict{String,String}())
end

@testset "Closure scatter data (headless)" begin
    ds = closure_scatter_fixture()
    programme = [0.5, -0.1, 0.0, 0.0]
    pairs = [(:BETA, :F3), (:BF, :F1)]
    data = closure_scatter_data(ds, pairs; programme)
    @test data.run_ids == ["BETA-F3", "BF-F1"]
    @test data.closure_pairs == [("BETA", "F3"), ("BF", "F1")]
    @test data.sectors == ["s1", "s2", "s3", "s4"]
    @test data.shocked == [true, true, false, false]
    @test data.output[:, 1] ≈ [6.25, -2.75, 0.85, -0.05]
    @test data.price[:, 1] ≈ [0.6, -0.3, 0.03, -0.06]
    @test data.price[:, 2] == zeros(4)
    @test data.domar_weights ≈ [0.04, 0.32, 0.01, 0.09]
    named = closure_scatter_data(ds,
        ["Alternative" => (:BETA, :F3), "Benchmark" => (:BF, :F1)]; programme)
    @test named.closure_labels == ["Alternative", "Benchmark"]
    @test named.closure_pairs == data.closure_pairs
    @test named.output == data.output && named.price == data.price
    mixed = closure_scatter_data(ds, ["Alternative" => (:BETA, :F3), (:BF, :F1)]; programme)
    @test mixed.closure_labels == ["Alternative", "BF–F1"]
    @test closure_scatter_data(ds, ["Named" => (:BF => :F1)]; programme).closure_labels == ["Named"]
    @test closure_scatter_data(ds, ["Named" => (:BF, :F1)]; programme,
        closure_labels = ["Override"]).closure_labels == ["Override"]
    for selection in ([" " => (:BF, :F1)], [123 => (:BF, :F1)],
            ["Same" => (:BF, :F1), "Same" => (:BETA, :F2)],
            ["First" => (:BF, :F1), "Second" => (:BF, :F1)])
        @test_throws ArgumentError closure_scatter_data(ds, selection; programme)
    end
    custom = closure_scatter_data(ds, [:BF => :F2, :BETA => :F1]; programme,
        domar_weights = [1, 4, 9, 16],
        sector_labels = ["A", "B", "C", "D"], closure_labels = ["First", "Second"])
    @test custom.labels == ["A", "B", "C", "D"]
    @test custom.closure_labels == ["First", "Second"]
    @test custom.domar_weights == [1, 4, 9, 16]
    @test size(custom.output) == size(custom.price) == (4, 2)
    @test !any(closure_scatter_data(ds, pairs; programme = zeros(4)).shocked)
    @test all(closure_scatter_data(ds, pairs; programme = ones(4)).shocked)
    @test_throws DimensionMismatch closure_scatter_data(ds, pairs; programme = [1.0])
    @test_throws DimensionMismatch closure_scatter_data(ds, pairs; programme, domar_weights = [1.0])
    for weights in ([1, -1, 2, 3], zeros(4), [NaN, 1, 2, 3], [Inf, 1, 2, 3])
        @test_throws ArgumentError closure_scatter_data(ds, pairs; programme, domar_weights = weights)
    end
    for value in (NaN, Inf)
        @test_throws ArgumentError closure_scatter_data(ds, pairs; programme = [value, 0, 0, 0])
        bad = deepcopy(ds); bad.cells[1].quantities[1] = value
        @test_throws ArgumentError closure_scatter_data(bad, pairs; programme)
        bad = deepcopy(ds); bad.cells[1].prices[1] = value
        @test_throws ArgumentError closure_scatter_data(bad, pairs; programme)
    end
    bad = deepcopy(ds); bad.cells[2].sectors[1:2] = ["s2", "s1"]
    @test_throws ArgumentError closure_scatter_data(bad, [(:BF, :F2)]; programme)
    bad = deepcopy(ds); pop!(bad.cells[1].prices)
    @test_throws DimensionMismatch closure_scatter_data(bad, pairs; programme)
    bad = deepcopy(ds); bad.cells[1].consumption .= NaN
    @test size(closure_scatter_data(bad, pairs; programme).output) == (4, 2)
    for selection in ([], [:BF], ["BF-F1"], [(:BF,)], [(1, 2)],
            [(:unknown, :F1)], [(:BF, :F1), (:BF, :F1)])
        @test_throws ArgumentError closure_scatter_data(ds, selection; programme)
    end
    ambiguous = deepcopy(ds); push!(ambiguous.cells, first(ambiguous.cells))
    @test_throws ArgumentError closure_scatter_data(ambiguous, [(:BF, :F1)]; programme)
    if Base.get_extension(OpenCGE, :OpenCGEGLMakieExt) === nothing
        @test_throws ErrorException plot_closure_scatter(ds, pairs; programme)
    end
end

if get(ENV, "OPENCGE_TEST_GLMAKIE", "false") == "true"
    using GLMakie
    GLMakie.activate!(; visible = false)

    @testset "Closure scatter graphics" begin
        ds = closure_scatter_fixture()
        programme = [0.5, -0.1, 0.0, 0.0]
        pairs = [(:BF, :F1), (:BETA, :F2), (:BETA, :F3)]
        weights = [1, 4, 9, 16]
        names = ["Benchmark", "Alternative", "Debt financed"]
        fig = plot_closure_scatter(ds, [name => pair for (name, pair) in zip(names, pairs)];
            programme, domar_weights = weights,
            closure_colors = [:red, :blue, :green])
        GLMakie.Makie.update_state_before_display!(fig)
        ax = only(filter(b -> b isa Axis, fig.content))
        @test size(ax.scene, 2) > 0.65 * size(fig.scene, 2)
        plots(f) = filter(p -> p isa GLMakie.Makie.Scatter,
            only(filter(b -> b isa Axis, f.content)).scene.plots)
        points = plots(fig)
        @test length(points) == 6
        @test all(length(p[1][]) == 2 for p in points)
        circle = GLMakie.Makie.to_spritemarker(:circle)
        @test all(p.marker[] == circle for p in points)
        data = closure_scatter_data(ds, pairs; programme)
        for k in 1:3
            @test first.(points[k][1][]) ≈ data.output[3:4, k]
            @test last.(points[k][1][]) ≈ data.price[3:4, k]
            @test first.(points[k + 3][1][]) ≈ data.output[1:2, k]
            @test last.(points[k + 3][1][]) ≈ data.price[1:2, k]
            # Areas, rather than diameters, carry the Domar-weight ratio.
            sizes = vcat(points[k + 3].markersize[], points[k].markersize[])
            @test all(s -> s[1] == s[2], sizes)
            diameters = first.(sizes)
            @test diameters ≈ [6, 12, 18, 24]
            @test diameters .^ 2 ./ diameters[1]^2 ≈ weights
            @test points[k].color[].alpha ≈ 0.35
            @test points[k + 3].color[].alpha == 1
            @test points[k].color[].r >= points[k + 3].color[].r
            @test points[k].color[].g >= points[k + 3].color[].g
            @test points[k].color[].b >= points[k + 3].color[].b
        end
        @test length(unique(p.color[] for p in points[4:6])) == 3
        legends = filter(b -> b isa Legend, fig.content)
        @test length(legends) == 2
        for legend in legends
            box = legend.layoutobservables.computedbbox[]
            @test all(box.origin .>= 0)
            @test all(box.origin .+ box.widths .<= size(fig.scene) .+ 1)
        end
        entries = only(first(legends).entrygroups[])[2]
        @test [entry.attributes[:label][] for entry in entries] == names
        for k in 1:3
            @test GLMakie.Makie.to_color(only(entries[k].elements).attributes[:markercolor][]) == points[k + 3].color[]
        end
        auto = plot_closure_scatter(ds, pairs; programme)
        permuted = plot_closure_scatter(ds, reverse(pairs); programme)
        @test [p.color[] for p in plots(auto)[1:3]] == reverse([p.color[] for p in plots(permuted)[1:3]])
        subset = plot_closure_scatter(ds, pairs[2:2]; programme)
        @test first(plots(subset)).color[] == plots(auto)[2].color[]
        shade = plot_closure_scatter(ds, pairs; programme = zeros(4))
        opaque = plot_closure_scatter(ds, pairs; programme = ones(4))
        @test first(plots(shade)).markersize[] == first(plots(opaque)).markersize[]
        @test length(plots(shade)) == length(plots(opaque)) == 3
        many = closure_scatter_fixture(; labours = ["BF", "ALPHA", "BETA", "GAMMA", "DELTA"])
        dense = plot_closure_scatter(many, [(c.labour, c.financing) for c in many.cells]; programme)
        @test length(unique(p.color[] for p in plots(dense)[1:15])) == 15
        @test all(p.marker[] == circle for p in plots(dense))
        denseax = only(filter(b -> b isa Axis, dense.content))
        @test denseax.xtickformat[] == "{:.1f}%"
        flat = plot_closure_scatter(ds, [(:BF, :F1)]; programme)
        flatax = only(filter(b -> b isa Axis, flat.content))
        GLMakie.Makie.update_state_before_display!(flat)
        @test flatax.finallimits[].widths[2] > 0
        @test all(all(iszero, last.(p[1][])) for p in plots(flat))
        @test_throws ArgumentError plot_closure_scatter(ds, pairs; programme, unshocked_alpha = 0)
        @test_throws ArgumentError plot_closure_scatter(ds, pairs; programme, max_markersize = -1)
        @test_throws DimensionMismatch plot_closure_scatter(ds, pairs; programme, closure_colors = [:red])
        @test_throws ArgumentError plot_closure_scatter(ds, pairs; programme, closure_colors = [:red, :red, :blue])
        mktempdir() do dir
            path = joinpath(dir, "scatter.png")
            save(path, fig)
            @test size(GLMakie.FileIO.load(path)) == (650, 1100)
        end
    end
end
