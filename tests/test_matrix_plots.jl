using OpenCGE
using Test

"""Synthetic sectoral responses for both manuscript-panel geometries."""
function matrix_panel_fixture()
	labels = ["Construction", "Electrical equipment", "Transport equipment",
		"Business services", "Basic metals", "Energy supply", "Chemicals",
		"Wholesale trade", "Food products", "Other services"]
	n = length(labels)
	y0 = collect(range(10.0, 5.0; length = n))
	baseline = MatrixBaseline(ones(n), y0, ones(n), y0 ./ 2, zeros(n), zeros(n))
	dq = [4.5, 3.8, 2.9, 0.8, 0.82, 0.84, -0.4, 0.2, 0.1, -0.3]
	dp = [0.8, 0.5, 0.6, 0.12, 0.11, 0.13, -0.08, 0.02, -0.01, -0.03]
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

	@testset "Print-ready matrix panels (GLMakie)" begin
		fixture = matrix_panel_fixture()
		(; ds, scatter, bars, variants, programme) = fixture
		for cell in (scatter, bars)
			fig = plot_matrix_panel(ds, cell.run_id; programme, variants, style = :manuscript)
			@test Tuple(size(fig.scene)) == (1100, 560)
			texts = [b.text[] for b in fig.content if b isa Label]
			@test all(t -> !any(s -> occursin(s, t),
				("Shock effect", "Source:", "5x3 matrix panel", cell.run_id)), texts)
			@test isempty(texts)
			legends = _mp_test_legends(fig)
			@test length(legends) == 2
			@test all(l -> l.labelsize[] >= 17, legends)
			for legend in legends
				bb = legend.layoutobservables.computedbbox[]
				@test all(bb.origin .>= 0)
				@test all(bb.origin .+ bb.widths .<= size(fig.scene))
			end
			left = only(filter(l -> any(t -> occursin("Programme demand", t),
				_mp_test_legend_labels(l)), legends))
			@test "Range over σ ∈ {0.6, 0.99}\n(including central cell)" in _mp_test_legend_labels(left)
			@test any(e -> e isa LineElement, last(left.entrygroups[][1][2]).elements)
			axes = filter(b -> b isa Axis, fig.content)
			@test all(a -> a.yticklabelsize[] >= 16 && a.xlabelsize[] >= 18 &&
				a.ylabelsize[] >= 18, axes)
			right = axes[2]
			texts = filter(p -> p isa GLMakie.Makie.Text, right.scene.plots)
			@test length(texts) == 6
			bounds = [GLMakie.Makie.boundingbox(p, :pixel) for p in texts]
			@test all(p -> p.fontsize[] >= 14, texts)
			@test all(bb -> all(bb.origin[1:2] .>= -1) &&
				all(bb.origin[1:2] .+ bb.widths[1:2] .<= size(right.scene) .+ 1), bounds)
			for i in eachindex(bounds), j in 1:(i - 1)
				a, b = bounds[i], bounds[j]
				overlap = min.(a.origin[1:2] .+ a.widths[1:2], b.origin[1:2] .+ b.widths[1:2]) .-
					max.(a.origin[1:2], b.origin[1:2])
				@test any(overlap .<= 0)
			end
		end

		fig = plot_matrix_panel(ds, bars.run_id; programme, style = :manuscript)
		labels = vcat(_mp_test_legend_labels.(_mp_test_legends(fig))...)
		@test "No σ variant converged" in labels
		@test !any(t -> occursin("Range over σ", t), labels)
		failed = Dict(bars.run_id => (low = nothing, high = nothing))
		fig = plot_matrix_panel(ds, bars.run_id; programme, variants = failed, style = :manuscript)
		@test "No σ variant converged" in vcat(_mp_test_legend_labels.(_mp_test_legends(fig))...)

		fig = plot_matrix_panel(ds, bars.run_id; programme)
		@test Tuple(size(fig.scene)) == (1900, 900)
		texts = [b.text[] for b in fig.content if b isa Label]
		@test any(t -> occursin("5x3 matrix panel", t), texts)
		@test any(t -> occursin("Source: " * bars.run_id, t), texts)
		fig = plot_matrix_panel(ds, bars.run_id; programme, title = "")
		@test !any(b -> b isa Label && occursin("5x3 matrix panel", b.text[]), fig.content)
		@test_throws ArgumentError plot_matrix_panel(ds, bars.run_id; programme, style = :unknown)
		mktempdir() do outdir
			paths = save_matrix_panels(ds; outdir, programme, variants, style = :manuscript)
			@test length(paths) == 2
			@test all(isfile, paths)
			@test all(p -> size(GLMakie.FileIO.load(p)) == (1120, 2200), paths)
			@test_throws ArgumentError save_matrix_panels(ds; outdir, programme, style = :unknown)
		end
	end
else
	@testset "Print-ready matrix panels (optional GLMakie)" begin
		@test_skip false
	end
end
