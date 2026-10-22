using OpenCGE
using Test

"""Synthetic sectoral responses for both manuscript-panel geometries."""
function matrix_panel_fixture(; extra_sectors::Int = 0)
	labels = ["Construction", "Electrical equipment", "Transport equipment",
		"Business services", "Basic metals", "Energy supply", "Chemicals",
		"Wholesale trade", "Food products", "Other services"]
	append!(labels, ["Other sector $j" for j in 1:extra_sectors])
	n = length(labels)
	y0 = collect(range(10.0, 5.0; length = n))
	baseline = MatrixBaseline(ones(n), y0, ones(n), y0 ./ 2, zeros(n), zeros(n))
	dq = [4.5, 3.8, 2.9, 0.8, 0.82, 0.84, -0.4, 0.2, 0.1, -0.3]
	dp = [0.8, 0.5, 0.6, 0.12, 0.11, 0.13, -0.08, 0.02, -0.01, -0.03]
	append!(dq, fill(0.05, extra_sectors))
	append!(dp, fill(0.01, extra_sectors))
	makecell(rid, labour, scale, moving) = MatrixCellData(rid, labour, "F1",
		1.0, 0.5, "ok", string.(1:n), labels,
		moving ? 1 .+ scale .* dp ./ 100 : ones(n),
		1 .+ scale .* dq ./ 100, ones(n), 1 .+ scale .* dq ./ 200,
		ones(n), ones(n), Dict{String,Any}(), Dict{String,Any}())
	scatter = makecell("synthetic-moving", "BETA", 1.0, true)
	bars = makecell("synthetic-flat", "ALPHA", 1.0, false)
	variants = Dict(c.run_id => (low = makecell(c.run_id * "-low", c.labour, 0.8,
		c === scatter), high = makecell(c.run_id * "-high", c.labour, 1.2,
		c === scatter)) for c in (scatter, bars))
	ds = MatrixDataset("synthetic-panels", baseline, [scatter, bars],
		Dict{String,Any}(), Dict{String,String}())
	programme = vcat([0.6, 0.4, 0.2], zeros(n - 3))
	return (; ds, scatter, bars, variants, programme)
end

# Opt in with OPENCGE_TEST_GLMAKIE=true and GLMakie on the load path.
# The ordinary kernel suite remains headless and synthetic-only.
if get(ENV, "OPENCGE_TEST_GLMAKIE", "false") == "true"
	using GLMakie
	GLMakie.activate!(; visible = false)

	_mp_test_legends(fig) = filter(b -> b isa Legend, fig.content)
	_mp_test_legend_labels(leg) = [entry.attributes[:label][]
		for (_, entries) in leg.entrygroups[] for entry in entries]

	_mp_test_axes(fig) = filter(b -> b isa Axis, fig.content)
	_mp_test_texts(ax) = filter(p -> p isa GLMakie.Makie.Text, ax.scene.plots)
	function _mp_test_nonoverlapping_text(ax)
		texts = _mp_test_texts(ax)
		bounds = [GLMakie.Makie.boundingbox(p, :pixel) for p in texts]
		@test all(p -> p.fontsize[] >= 17, texts)
		@test all(bb -> all(bb.origin[1:2] .>= -1) &&
			all(bb.origin[1:2] .+ bb.widths[1:2] .<= size(ax.scene) .+ 1), bounds)
		for i in eachindex(bounds), j in 1:(i - 1)
			a, b = bounds[i], bounds[j]
			overlap = min.(a.origin[1:2] .+ a.widths[1:2], b.origin[1:2] .+ b.widths[1:2]) .-
				max.(a.origin[1:2], b.origin[1:2])
			@test any(overlap .<= 0)
		end
	end

	@testset "Print-ready matrix panels (GLMakie)" begin
		(; ds, scatter, bars, variants, programme) = matrix_panel_fixture()
		curated = ["Short $j" for j in eachindex(programme)]
		for cell in (scatter, bars)
			fig = plot_matrix_panel(ds, cell.run_id; programme, variants,
				sector_labels = curated, style = :manuscript, title = "Do not print this title")
			@test Tuple(size(fig.scene)) == (1100, 520)
			texts = [b.text[] for b in fig.content if b isa Label]
			@test all(t -> !any(s -> occursin(s, t),
				("Shock effect", "Source:", "5x3 matrix panel", cell.run_id, "Do not print this title")), texts)
			@test isempty(texts)
			legend = only(_mp_test_legends(fig))
			@test legend.labelsize[] >= 17
			@test legend.orientation[] == :horizontal
			@test legend.nbanks[] <= 2
			bb = legend.layoutobservables.computedbbox[]
			@test all(bb.origin .>= 0)
			@test all(bb.origin .+ bb.widths .<= size(fig.scene))
			labels = _mp_test_legend_labels(legend)
			@test "σ ∈ {0.6, 0.99} range" in labels
			@test ("Price" in labels) == (cell === scatter)
			@test any(e -> e isa LineElement,
				only(filter(e -> e.attributes[:label][] == "σ ∈ {0.6, 0.99} range",
					legend.entrygroups[][1][2])).elements)
			axes = _mp_test_axes(fig)
			@test all(a -> a.xticklabelsize[] >= 16 && a.yticklabelsize[] >= 16 &&
				a.xlabelsize[] >= 17 && a.ylabelsize[] >= 17, axes)
			@test all(a -> size(a.scene, 2) >= 0.70 * size(fig.scene, 2), axes)
			@test axes[1].yticks[][2] == ["$k $(curated[j])" for (k, j) in enumerate(1:3)]
			@test axes[1].yreversed[]
			@test axes[1].xticklabelrotation[] == 0
			@test axes[1].xtickformat[] == "{:.0f}%"
			@test axes[2].ytickformat[] == (cell === scatter ? "{:.1f}%" : "{:.0f}%")
			_mp_test_nonoverlapping_text(axes[2])
			annotations = [only(p.text[]) for p in _mp_test_texts(axes[2])]
			@test all(string(k) in annotations for k in 1:3)
			@test count(t -> !all(isdigit, t), annotations) <= 2
			@test_throws DimensionMismatch plot_matrix_panel(ds, cell.run_id;
				programme, style = :manuscript, sector_labels = curated[1:2])
		end

		# All seven programme sectors are numbered, including a tight scatter cluster.
		seven = vcat(collect(0.7:-0.1:0.1), zeros(3))
		for cell in (scatter, bars)
			fig = plot_matrix_panel(ds, cell.run_id; programme = seven, style = :manuscript)
			ax = _mp_test_axes(fig)[2]
			@test all(string(k) in [only(p.text[]) for p in _mp_test_texts(ax)] for k in 1:7)
			_mp_test_nonoverlapping_text(ax)
		end

		# Narrow bars at 71-sector print density still show every programme number.
		dense = matrix_panel_fixture(; extra_sectors = 61)
		fig = plot_matrix_panel(dense.ds, dense.bars.run_id;
			programme = vcat(seven, zeros(61)), style = :manuscript)
		ax = _mp_test_axes(fig)[2]
		@test all(string(k) in [only(p.text[]) for p in _mp_test_texts(ax)] for k in 1:7)
		_mp_test_nonoverlapping_text(ax)

		ext = Base.get_extension(OpenCGE, :OpenCGEGLMakieExt)
		@test ext._mp_panel_label("Wood.cork_exc.furniture.and.other.products") == "Wood cork exc furniture…"
		@test ext._mp_panel_label("Glass") == "Glass"
		@test ext._mp_panel_tickformat(0.0, 2.99) == "{:.1f}%"
		@test ext._mp_panel_tickformat(0.0, 3.0) == "{:.0f}%"
		# Dense markers leave no safe outlier label box: drop it.
		fig = Figure(size = (500, 500))
		ax = Axis(fig[1, 1]); xlims!(ax, 0, 1); ylims!(ax, 0, 1)
		GLMakie.Makie.update_state_before_display!(fig)
		points = [(x, y) for x in 0:0.01:1 for y in 0:0.01:1]
		@test !ext._mp_panel_annotation!(ax, "Cannot fit", 0.5, 0.5,
			(0., 1., 0., 1.), points, NTuple{4,Float64}[])
		@test isempty(_mp_test_texts(ax))

		for vars in (Dict{String,Any}(), Dict(bars.run_id => (low = nothing, high = nothing)))
			fig = plot_matrix_panel(ds, bars.run_id; programme, variants = vars, style = :manuscript)
			labels = _mp_test_legend_labels(only(_mp_test_legends(fig)))
			@test !("Price" in labels)
			@test !any(t -> occursin("range", t), labels)
		end
		fig = plot_matrix_panel(ds, bars.run_id; programme, right = :scatter, style = :manuscript)
		@test _mp_test_axes(fig)[2].subtitle[] == "Prices at baseline"

		fig = plot_matrix_panel(ds, bars.run_id; programme)
		@test Tuple(size(fig.scene)) == (1900, 900)
		texts = [b.text[] for b in fig.content if b isa Label]
		@test any(t -> occursin("5x3 matrix panel", t), texts)
		@test any(t -> occursin("Source: " * bars.run_id, t), texts)
		@test length(_mp_test_legends(fig)) == 2
		@test _mp_test_axes(fig)[1].xticks[][2] == bars.labels[1:3]
		fig = plot_matrix_panel(ds, bars.run_id; programme, title = "")
		@test !any(b -> b isa Label && occursin("5x3 matrix panel", b.text[]), fig.content)
		@test_throws ArgumentError plot_matrix_panel(ds, bars.run_id; programme, style = :unknown)
		mktempdir() do outdir
			paths = save_matrix_panels(ds; outdir, programme, variants,
				sector_labels = curated, style = :manuscript)
			@test length(paths) == 2
			@test all(isfile, paths)
			@test all(p -> size(GLMakie.FileIO.load(p)) == (1040, 2200), paths)
			@test_throws DimensionMismatch save_matrix_panels(ds; outdir, programme,
				style = :manuscript, sector_labels = curated[1:2])
			@test_throws ArgumentError save_matrix_panels(ds; outdir, programme, style = :unknown)
		end
	end
else
	@testset "Print-ready matrix panels (optional GLMakie)" begin
		@test_skip false
	end
end
