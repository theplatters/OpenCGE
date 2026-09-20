using Printf

# 5x3 evaluation-matrix figures (WP2; exploratory + comparable across closures).
#
# All panels use the canonical cell ordering
# (`matrix_labour_order()` x `matrix_financing_order()`), percent units with
# labelled axes, and one shared colour scale per metric per figure. Sectoral
# panels sort sectors by baseline gross output descending. Every figure states
# its baseline (the design reference / no-programme baseline) in a footnote
# label. Missing cells render as NaN (blank); NaN sectoral relatives (zero
# baselines) are masked everywhere. All functions return a `Figure`;
# `save_matrix_figures` writes files and returns the paths.

"""Footnote naming the baseline every matrix figure is measured against."""
const _MP_BASELINE_NOTE =
	"Shock effect relative to the design reference (no-programme baseline); 1.0 = unchanged."

"""Canonical-order sort key for a dataset cell (unknown ids sort last)."""
function _mp_cell_key(cell::BeyondHulten.MatrixCellData)
	labour = BeyondHulten.matrix_labour_order()
	financing = BeyondHulten.matrix_financing_order()
	li = findfirst(==(cell.labour), labour)
	fi = findfirst(==(cell.financing), financing)
	return (li === nothing ? typemax(Int) : li, fi === nothing ? typemax(Int) : fi)
end

"""Cells of `ds` in canonical labour x financing order."""
function _mp_sorted_cells(ds::BeyondHulten.MatrixDataset)
	return sort(collect(ds.cells); by = _mp_cell_key)
end

"""Short `labour-financing` label for a cell."""
_mp_cell_label(cell::BeyondHulten.MatrixCellData) = cell.labour * "-" * cell.financing

"""Aggregate metric as `Float64` (`NaN` when missing or non-numeric)."""
function _mp_metric(cell::BeyondHulten.MatrixCellData, key::AbstractString)
	v = get(cell.metrics, String(key), NaN)
	v isa Real || return NaN
	return Float64(v)
end

"""Diagnostic entry as `Float64` (`NaN` when missing or non-numeric)."""
function _mp_diag(cell::BeyondHulten.MatrixCellData, key::AbstractString)
	v = get(cell.diagnostics, String(key), NaN)
	v isa Real || return NaN
	return Float64(v)
end

"""Percent annotation for a heatmap cell (`"n/a"` for missing values)."""
function _mp_pct_text(x::Real)
	isfinite(x) || return "n/a"
	abs(Float64(x)) < 0.005 && return "0.00%"
	return @sprintf("%.2f%%", x)
end

"""Figure title: explicit `title`, else a default naming the design."""
function _mp_title(title, default::AbstractString)
	title === nothing && return default
	return String(title)
end

"""Symmetric colour range around 0 over the finite entries (fallback ±1)."""
function _mp_symmetric_range(vals::AbstractVector{<:Real})
	f = filter(isfinite, vals)
	isempty(f) && return (-1.0, 1.0)
	m = maximum(abs, f)
	(!isfinite(m) || m == 0.0) && return (-1.0, 1.0)
	return (-m, m)
end

"""Quantile of a sorted vector at `q` (linear interpolation)."""
function _mp_quantile_sorted(s::AbstractVector{<:Real}, q::Real)
	n = length(s)
	n == 1 && return Float64(s[1])
	h = (n - 1) * clamp(Float64(q), 0.0, 1.0) + 1.0
	lo = clamp(floor(Int, h), 1, n)
	hi = clamp(ceil(Int, h), 1, n)
	lo == hi && return Float64(s[lo])
	t = h - lo
	return (1 - t) * Float64(s[lo]) + t * Float64(s[hi])
end

"""Shared symmetric robust range: `±max(|p2|, |p98|)` over the finite entries.

Falls back to the max-abs range when the percentile half-width is 0 or
non-finite. Heatmaps use this so single small-baseline outliers do not
dominate the colour scale; the box/strip summary keeps the full range."""
function _mp_robust_range(vals::AbstractVector{<:Real})
	f = filter(isfinite, vals)
	isempty(f) && return (-1.0, 1.0)
	s = sort(Float64.(f))
	p2 = _mp_quantile_sorted(s, 0.02)
	p98 = _mp_quantile_sorted(s, 0.98)
	m = max(abs(p2), abs(p98))
	(!isfinite(m) || m == 0.0) && return _mp_symmetric_range(vals)
	return (-m, m)
end

"""Data-limit colour range over the finite entries (fallback ±1)."""
function _mp_data_range(vals::AbstractVector{<:Real})
	f = filter(isfinite, vals)
	isempty(f) && return (-1.0, 1.0)
	lo, hi = minimum(f), maximum(f)
	(!isfinite(lo) || !isfinite(hi) || lo == hi) && return (lo - 1.0, hi + 1.0)
	return (lo, hi)
end

"""Sector positions sorting baseline gross output descending."""
function _mp_sector_order(ds::BeyondHulten.MatrixDataset)
	go = ds.baseline.prices .* ds.baseline.quantities
	return sortperm(go; rev = true)
end

"""Matrix of sectoral percent changes: rows = cells, cols = sector-order positions."""
function _mp_sector_matrix(cells::Vector{BeyondHulten.MatrixCellData},
		order::Vector{Int}, getvec::Function)
	nc = length(cells)
	ns = length(order)
	mat = fill(NaN, nc, ns)
	for (k, cell) in enumerate(cells)
		v = getvec(cell)
		length(v) == length(order) || continue
		for (sp, j) in enumerate(order)
			r = j <= length(v) ? v[j] : NaN
			mat[k, sp] = isfinite(r) ? 100.0 * (r - 1.0) : NaN
		end
	end
	return mat
end

"""Finite per-cell values for a box/strip summary (one vector per cell)."""
function _mp_box_values(mat::Matrix{Float64})
	return [filter(isfinite, mat[k, :]) for k in axes(mat, 1)]
end

"""Jittered x positions around integer group `k` (deterministic)."""
function _mp_jitter(k::Int, n::Int; width::Float64 = 0.3)
	n == 1 && return [Float64(k)]
	return [k + width * (2.0 * ((i * 0.6180339887498949) % 1.0) - 1.0) for i in 1:n]
end

"""Edge-aware label placement for scatter panels.

Returns `(align, offset)` for a point label: points in the right ~15% of
the x-range are labelled to the LEFT of the marker (`(:right, :bottom)`
with a negative x offset, so the text extends inward), mirrored for
points in the top ~15% of the y-range (label BELOW the marker); anywhere
else the existing alternating offsets are kept (`k` = cell position)."""
function _mp_scatter_align_offset(x::Real, y::Real,
		xlo::Real, xhi::Real, ylo::Real, yhi::Real, k::Int)
	xr = xhi - xlo
	yr = yhi - ylo
	x_edge = isfinite(xr) && xr > 0.0 && x >= xhi - 0.15 * xr
	y_edge = isfinite(yr) && yr > 0.0 && y >= yhi - 0.15 * yr
	if x_edge && y_edge
		return ((:right, :top), (-6, -6))
	elseif x_edge
		return ((:right, :bottom), (-6, 6))
	elseif y_edge
		return ((:left, :top), (6, -6))
	end
	if iseven(k)
		return ((:left, :bottom), (6, 6))
	else
		return ((:right, :top), (-6, -6))
	end
end

"""Expand `ax` limits ~7% beyond the data range so point labels fit."""
function _mp_pad_axis!(ax, xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}; pad::Float64 = 0.07)
	(isempty(xs) || isempty(ys)) && return nothing
	xlo, xhi = minimum(xs), maximum(xs)
	ylo, yhi = minimum(ys), maximum(ys)
	xr = xhi - xlo
	yr = yhi - ylo
	xr = xr == 0.0 ? max(abs(xlo), 1.0) : xr
	yr = yr == 0.0 ? max(abs(ylo), 1.0) : yr
	xlims!(ax, xlo - pad * xr, xhi + pad * xr)
	ylims!(ax, ylo - pad * yr, yhi + pad * yr)
	return nothing
end

function BeyondHulten.plot_matrix_overview(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1150))
	labs = BeyondHulten.matrix_labour_order()
	fins = BeyondHulten.matrix_financing_order()
	bykey = Dict{Tuple{String,String},BeyondHulten.MatrixCellData}(
		(c.labour, c.financing) => c for c in ds.cells)
	panels = [
		("Consumption, rel. change (%)", c -> 100.0 * _mp_metric(c, "consumption_rel"), true),
		("Employment, change (%)", c -> 100.0 * (_mp_metric(c, "employment") - 1.0), true),
		("Real GDP, rel. change (%)", c -> 100.0 * _mp_metric(c, "gdp_rel"), true),
		("GDP deflator, change (%)", c -> 100.0 * (_mp_metric(c, "gdp_deflator") - 1.0), true),
		("External position (% of GDP)", c -> 100.0 * _mp_metric(c, "external_position"), false),
		("Wage, change (%)", c -> 100.0 * (_mp_metric(c, "wage") - 1.0), true),
	]
	fig = Figure(size = size)
	Label(fig[1, 1:3], _mp_title(title, "5x3 matrix overview ($(ds.design))"), fontsize = 20)
	for (k, (pt, fn, centred)) in enumerate(panels)
		i = div(k - 1, 3) + 2
		j = mod(k - 1, 3) + 1
		g = GridLayout(fig[i, j])
		ax = Axis(g[1, 1]; title = pt, titlesize = 14,
			xticks = (1:length(fins), fins), yticks = (1:length(labs), labs),
			yreversed = true, xticklabelsize = 11, yticklabelsize = 11)
		M = fill(NaN, length(fins), length(labs))
		for (li, l) in enumerate(labs), (fi, f) in enumerate(fins)
			haskey(bykey, (l, f)) || continue
			M[fi, li] = fn(bykey[(l, f)])
		end
		cr = centred ? _mp_symmetric_range(vec(M)) : _mp_data_range(vec(M))
		hm = heatmap!(ax, 1:length(fins), 1:length(labs), M;
			colormap = :RdBu, colorrange = cr)
		m = maximum((abs(cr[1]), abs(cr[2])))
		for (li, l) in enumerate(labs), (fi, f) in enumerate(fins)
			v = M[fi, li]
			isfinite(v) || continue
			tc = (centred && abs(v) > 0.55 * m) ? :white : :black
			text!(ax, fi, li; text = _mp_pct_text(v),
				align = (:center, :center), fontsize = 11, color = tc)
		end
		Colorbar(g[1, 2], hm; label = "%")
	end
	Label(fig[4, 1:3], _MP_BASELINE_NOTE, fontsize = 13)
	return fig
end

function BeyondHulten.plot_matrix_wages(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1150))
	labs = BeyondHulten.matrix_labour_order()
	fins = BeyondHulten.matrix_financing_order()
	li_of(c) = something(findfirst(==(c.labour), labs), 0)
	fi_of(c) = something(findfirst(==(c.financing), fins), 0)
	cells = _mp_sorted_cells(ds)
	colors = Makie.wong_colors()
	fig = Figure(size = size)
	Label(fig[1, 1:2], _mp_title(title, "5x3 matrix wages ($(ds.design))"), fontsize = 20)
	ax1 = Axis(fig[2, 1:2]; title = "(a) Aggregate wage change by cell (wage-bill-weighted for sectoral-wage cells)",
		ylabel = "change (%)", ytickformat = "{:.2f}%",
		xticks = (1:length(labs), labs))
	xs = Int[]
	dg = Int[]
	hs = Float64[]
	cs = []
	seen_fin = Int[]
	for c in cells
		(li, fi) = (li_of(c), fi_of(c))
		(li == 0 || fi == 0) && continue
		v = 100.0 * (_mp_metric(c, "wage") - 1.0)
		isfinite(v) || continue
		push!(xs, li)
		push!(dg, fi)
		push!(hs, v)
		push!(cs, colors[fi])
		fi in seen_fin || push!(seen_fin, fi)
	end
	if !isempty(hs)
		barplot!(ax1, xs, hs; dodge = dg, color = cs)
		sort!(seen_fin)
		elements = [PolyElement(polycolor = colors[fi]) for fi in seen_fin]
		axislegend(ax1, elements, [fins[fi] for fi in seen_fin];
			position = :rt, labelsize = 12)
	else
		text!(ax1, 0.5, 0.5; text = "No wage metric available.",
			align = (:center, :center), fontsize = 13)
	end
	ax2 = Axis(fig[3, 1]; title = "(b) Sectoral wage changes (non-degenerate cells)",
		ylabel = "change (%)", ytickformat = "{:.2f}%")
	dispersion(c) = begin
		f = filter(isfinite, c.wages)
		isempty(f) ? 0.0 : maximum(f) - minimum(f)
	end
	qual = [c for c in cells if dispersion(c) > 1e-12]
	if isempty(qual)
		xlims!(ax2, 0, 1)
		ylims!(ax2, 0, 1)
		text!(ax2, 0.5, 0.5;
			text = "All cells carry a single (flat) wage — sectoral panel skipped.",
			align = (:center, :center), fontsize = 13)
	else
		for (k, c) in enumerate(qual)
			fi_of(c) == 0 && continue
			vals = [100.0 * (w - 1.0) for w in c.wages if isfinite(w)]
			isempty(vals) && continue
			scatter!(ax2, _mp_jitter(k, length(vals)), vals;
				color = colors[fi_of(c)], markersize = 8)
		end
		ax2.xticks = (1:length(qual), [_mp_cell_label(c) for c in qual])
		ax2.xticklabelrotation = π / 4
	end
	ax3 = Axis(fig[3, 2]; title = "(c) Wage vs employment",
		xlabel = "wage change (%)", ylabel = "employment change (%)",
		xtickformat = "{:.2f}%", ytickformat = "{:.2f}%")
	hlines!(ax3, [0.0]; color = :gray, linestyle = :dash)
	vlines!(ax3, [0.0]; color = :gray, linestyle = :dash)
	scx = Float64[]
	scy = Float64[]
	pts = Tuple{Int,Float64,Float64,String,Int}[]
	for (k, c) in enumerate(cells)
		fi_of(c) == 0 && continue
		x = 100.0 * (_mp_metric(c, "wage") - 1.0)
		y = 100.0 * (_mp_metric(c, "employment") - 1.0)
		(isfinite(x) && isfinite(y)) || continue
		push!(pts, (k, x, y, _mp_cell_label(c), fi_of(c)))
		push!(scx, x)
		push!(scy, y)
	end
	# Coincident cells (e.g. ALPHA-F2 ≡ BETA-F2) share one marker and a
	# merged "A = B" label so the labels do not overprint; the identity
	# itself is a finding and stays visible.
	wgroups = Dict{Tuple{Float64,Float64},Vector{Tuple{Int,Float64,Float64,String,Int}}}()
	worder = Tuple{Float64,Float64}[]
	for p in pts
		key = (round(p[2]; digits = 6), round(p[3]; digits = 6))
		if !haskey(wgroups, key)
			wgroups[key] = Tuple{Int,Float64,Float64,String,Int}[]
			push!(worder, key)
		end
		push!(wgroups[key], p)
	end
	for key in worder
		members = wgroups[key]
		(k, x, y, _, fi) = members[1]
		lab = join([m[4] for m in members], " = ")
		scatter!(ax3, [x], [y]; color = colors[fi], markersize = 12)
		al, off = _mp_scatter_align_offset(x, y,
			minimum(scx), maximum(scx), minimum(scy), maximum(scy), k)
		text!(ax3, x, y; text = lab, align = al, fontsize = 10, offset = off)
	end
	if !isempty(scx)
		_mp_pad_axis!(ax3, scx, scy)
	end
	Label(fig[4, 1:2],
		_MP_BASELINE_NOTE * "  The CPI numeraire is 1 at every solution, so the recorded wages are real wages; panel (a) aggregates the sectoral-wage cells by wage bill." *
			"  Dispersion threshold for panel (b): 1e-12.",
		fontsize = 13)
	return fig
end

"""Shared sectors-x-cells heatmap figure backing prices and quantities."""
function _mp_sector_figure(ds::BeyondHulten.MatrixDataset, getvec::Function,
		default_title::AbstractString, panel_a::AbstractString, panel_b::AbstractString;
		title = nothing, size = (1800, 1500), top_labels::Int = 10, robust::Bool = true)
	cells = _mp_sorted_cells(ds)
	nc = length(cells)
	ns = length(ds.baseline.prices)
	order = _mp_sector_order(ds)
	mat = _mp_sector_matrix(cells, order, getvec)
	cr = robust ? _mp_robust_range(vec(mat)) : _mp_symmetric_range(vec(mat))
	clip_note = robust ?
		"shared scale, clipped at p98 = " * @sprintf("%.2f%%", cr[2]) :
		"shared scale, max |change| = " * @sprintf("%.2f%%", cr[2])
	fig = Figure(size = size)
	Label(fig[1, 1:2], _mp_title(title, default_title), fontsize = 20)
	g = GridLayout(fig[2, 1:2])
	rowsize!(fig.layout, 2, Makie.Fixed(900))
	labels = [_mp_cell_label(c) for c in cells]
	top = min(top_labels, ns)
	yticklabels = nc == 0 ? string.(1:top) : [cells[1].labels[order[k]] for k in 1:top]
	ax = Axis(g[1, 1]; title = panel_a * " (" * clip_note * ")",
		xticks = nc == 0 ? ([1], [""]) : (1:nc, labels),
		yticks = (1:top, yticklabels),
		yreversed = true, xticklabelrotation = π / 4,
		xticklabelsize = 11, yticklabelsize = 9)
	# heatmap(x, y, z) needs size(z) == (length(x), length(y)): mat[k, sp]
	# is already (cell, sector-position).
	hm = heatmap!(ax, 1:max(nc, 1), 1:ns,
		nc == 0 ? fill(NaN, 1, ns) : mat;
		colormap = :RdBu, colorrange = cr)
	Colorbar(g[1, 2], hm; label = "%")
	ax2 = Axis(fig[3, 1:2]; title = panel_b,
		ylabel = "change (%)", ytickformat = "{:.2f}%",
		xticks = nc == 0 ? ([1], [""]) : (1:nc, labels),
		xticklabelrotation = π / 4)
	bxs = Int[]
	bys = Float64[]
	for (k, vals) in enumerate(_mp_box_values(mat))
		isempty(vals) && continue
		append!(bxs, fill(k, length(vals)))
		append!(bys, vals)
	end
	if !isempty(bys)
		boxplot!(ax2, bxs, bys)
		scatter!(ax2, [x + 0.12 * (2.0 * ((i * 0.6180339887498949) % 1.0) - 1.0)
			for (i, x) in enumerate(bxs)], bys;
			markersize = 6, color = (:black, 0.35))
		hlines!(ax2, [0.0]; color = :gray, linestyle = :dash)
	else
		text!(ax2, 0.5, 0.5; text = "No finite sectoral values.",
			align = (:center, :center), fontsize = 13)
	end
	Label(fig[4, 1:2],
		_MP_BASELINE_NOTE * "  Sectors sorted by baseline gross output, descending." *
			"  Heatmap " * clip_note * "; box panel shows the full range.",
		fontsize = 13)
	return fig
end

function BeyondHulten.plot_matrix_prices(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1500), top_labels::Int = 10)
	return _mp_sector_figure(ds, c -> c.prices,
		"5x3 matrix prices ($(ds.design))",
		"(a) Sectoral price changes, %",
		"(b) Per-cell distribution of sectoral price changes, %";
		title = title, size = size, top_labels = top_labels, robust = false)
end

function BeyondHulten.plot_matrix_quantities(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1500), top_labels::Int = 10)
	return _mp_sector_figure(ds, c -> c.quantities,
		"5x3 matrix quantities ($(ds.design))",
		"(a) Sectoral quantity changes, %",
		"(b) Per-cell distribution of sectoral quantity changes, %";
		title = title, size = size, top_labels = top_labels, robust = false)
end

function BeyondHulten.plot_matrix_consumption(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1500), top_labels::Int = 10)
	labs = BeyondHulten.matrix_labour_order()
	fins = BeyondHulten.matrix_financing_order()
	li_of(c) = something(findfirst(==(c.labour), labs), 0)
	fi_of(c) = something(findfirst(==(c.financing), fins), 0)
	cells = _mp_sorted_cells(ds)
	nc = length(cells)
	ns = length(ds.baseline.prices)
	colors = Makie.wong_colors()
	fig = Figure(size = size)
	Label(fig[1, 1:2], _mp_title(title, "5x3 matrix consumption ($(ds.design))"), fontsize = 20)
	ax1 = Axis(fig[2, 1:2]; title = "(a) Aggregate consumption change by cell",
		ylabel = "change (%)", ytickformat = "{:.2f}%",
		xticks = (1:length(labs), labs))
	xs = Int[]
	dg = Int[]
	hs = Float64[]
	cs = []
	seen_fin = Int[]
	for c in cells
		(li, fi) = (li_of(c), fi_of(c))
		(li == 0 || fi == 0) && continue
		v = 100.0 * _mp_metric(c, "consumption_rel")
		isfinite(v) || continue
		push!(xs, li)
		push!(dg, fi)
		push!(hs, v)
		push!(cs, colors[fi])
		fi in seen_fin || push!(seen_fin, fi)
	end
	if !isempty(hs)
		barplot!(ax1, xs, hs; dodge = dg, color = cs)
		sort!(seen_fin)
		elements = [PolyElement(polycolor = colors[fi]) for fi in seen_fin]
		axislegend(ax1, elements, [fins[fi] for fi in seen_fin];
			position = :rt, labelsize = 12)
	else
		text!(ax1, 0.5, 0.5; text = "No consumption metric available.",
			align = (:center, :center), fontsize = 13)
	end
	order = _mp_sector_order(ds)
	mat = _mp_sector_matrix(cells, order, c -> c.consumption)
	cr = _mp_robust_range(vec(mat))
	clip_note_c = "shared scale, clipped at p98 = ±" * @sprintf("%.2f%%", cr[2])
	g = GridLayout(fig[3, 1])
	rowsize!(fig.layout, 3, Makie.Fixed(900))
	labels = [_mp_cell_label(c) for c in cells]
	top = min(top_labels, ns)
	yticklabels = nc == 0 ? string.(1:top) : [cells[1].labels[order[k]] for k in 1:top]
	ax2 = Axis(g[1, 1]; title = "(b) Sectoral consumption changes, % (" * clip_note_c * ")",
		xticks = nc == 0 ? ([1], [""]) : (1:nc, labels),
		yticks = (1:top, yticklabels),
		yreversed = true, xticklabelrotation = π / 4,
		xticklabelsize = 10, yticklabelsize = 9, titlesize = 14)
	hm = heatmap!(ax2, 1:max(nc, 1), 1:ns,
		nc == 0 ? fill(NaN, 1, ns) : mat;
		colormap = :RdBu, colorrange = cr)
	Colorbar(g[1, 2], hm; label = "%")
	ax3 = Axis(fig[3, 2]; title = "(c) Consumption vs employment",
		xlabel = "consumption change (%)", ylabel = "employment change (%)",
		xtickformat = "{:.2f}%", ytickformat = "{:.2f}%")
	hlines!(ax3, [0.0]; color = :gray, linestyle = :dash)
	vlines!(ax3, [0.0]; color = :gray, linestyle = :dash)
	sccx = Float64[]
	sccy = Float64[]
	cpts = Tuple{Int,Float64,Float64,String,Int}[]
	for (k, c) in enumerate(cells)
		fi_of(c) == 0 && continue
		x = 100.0 * _mp_metric(c, "consumption_rel")
		y = 100.0 * (_mp_metric(c, "employment") - 1.0)
		(isfinite(x) && isfinite(y)) || continue
		push!(cpts, (k, x, y, _mp_cell_label(c), fi_of(c)))
		push!(sccx, x)
		push!(sccy, y)
	end
	# Coincident cells share one marker and a merged "A = B" label so the
	# labels do not overprint; the identity itself is a finding and stays
	# visible.
	cgroups = Dict{Tuple{Float64,Float64},Vector{Tuple{Int,Float64,Float64,String,Int}}}()
	corder = Tuple{Float64,Float64}[]
	for p in cpts
		key = (round(p[2]; digits = 6), round(p[3]; digits = 6))
		if !haskey(cgroups, key)
			cgroups[key] = Tuple{Int,Float64,Float64,String,Int}[]
			push!(corder, key)
		end
		push!(cgroups[key], p)
	end
	for key in corder
		members = cgroups[key]
		(k, x, y, _, fi) = members[1]
		lab = join([m[4] for m in members], " = ")
		scatter!(ax3, [x], [y]; color = colors[fi], markersize = 12)
		al, off = _mp_scatter_align_offset(x, y,
			minimum(sccx), maximum(sccx), minimum(sccy), maximum(sccy), k)
		text!(ax3, x, y; text = lab, align = al, fontsize = 10, offset = off)
	end
	if !isempty(sccx)
		_mp_pad_axis!(ax3, sccx, sccy)
	end
	Label(fig[4, 1:2], _MP_BASELINE_NOTE * "  Heatmap " * clip_note_c * "; scatter panel (c) keeps the full range.", fontsize = 13)
	return fig
end

function BeyondHulten.plot_matrix_trade(ds::BeyondHulten.MatrixDataset;
		title = nothing, size = (1800, 1500), top_labels::Int = 10)
	cells = _mp_sorted_cells(ds)
	nc = length(cells)
	ns = length(ds.baseline.prices)
	colors = Makie.wong_colors()
	fig = Figure(size = size)
	Label(fig[1, 1:2], _mp_title(title, "5x3 matrix trade ($(ds.design))"), fontsize = 20)
	gross_series = ["Exports X", "Total imports M"]
	net_series = ["Net external position", "Programme financing"]
	cellticks = nc == 0 ? ([1], [""]) : (1:nc, [_mp_cell_label(c) for c in cells])
	ax1a = Axis(fig[2, 1]; title = "(a1) Gross flows: exports X and total imports M",
		ylabel = "% of GDP", ytickformat = "{:.2f}%",
		xticks = cellticks,
		xticklabelrotation = π / 4)
	ax1b = Axis(fig[2, 2]; title = "(a2) Net external position and programme financing",
		ylabel = "% of GDP", ytickformat = "{:.2f}%",
		xticks = cellticks,
		xticklabelrotation = π / 4)
	gx = Int[]
	gd = Int[]
	gh = Float64[]
	gc = []
	nx = Int[]
	nd = Int[]
	nh = Float64[]
	nc2 = []
	for (k, c) in enumerate(cells)
		gross = [
			100.0 * _mp_diag(c, "gdp_x"),
			-100.0 * (_mp_diag(c, "gdp_m_final") + _mp_diag(c, "gdp_m_int")),
		]
		for (si, v) in enumerate(gross)
			isfinite(v) || continue
			push!(gx, k)
			push!(gd, si)
			push!(gh, v)
			push!(gc, colors[si])
		end
		net = [
			100.0 * _mp_metric(c, "external_position"),
			100.0 * _mp_metric(c, "programme_financing"),
		]
		for (si, v) in enumerate(net)
			isfinite(v) || continue
			push!(nx, k)
			push!(nd, si)
			push!(nh, v)
			push!(nc2, colors[si + 2])
		end
	end
	if !isempty(gh)
		barplot!(ax1a, gx, gh; dodge = gd, color = gc)
		elements = [PolyElement(polycolor = colors[si]) for si in 1:length(gross_series)]
		axislegend(ax1a, elements, gross_series; position = :rt, labelsize = 12)
		hlines!(ax1a, [0.0]; color = :gray, linestyle = :dash)
	else
		text!(ax1a, 0.5, 0.5; text = "No gross-flow decomposition available.",
			align = (:center, :center), fontsize = 13)
	end
	if !isempty(nh)
		barplot!(ax1b, nx, nh; dodge = nd, color = nc2)
		elements = [PolyElement(polycolor = colors[si + 2]) for si in 1:length(net_series)]
		axislegend(ax1b, elements, net_series; position = :rt, labelsize = 12)
		hlines!(ax1b, [0.0]; color = :gray, linestyle = :dash)
	else
		text!(ax1b, 0.5, 0.5; text = "No net-position decomposition available.",
			align = (:center, :center), fontsize = 13)
	end
	order = _mp_sector_order(ds)
	mat = _mp_sector_matrix(cells, order, c -> c.imports)
	cr = _mp_robust_range(vec(mat))
	clip_note_t = "shared scale, clipped at p98 = ±" * @sprintf("%.2f%%", cr[2])
	g = GridLayout(fig[3, 1:2])
	rowsize!(fig.layout, 3, Makie.Fixed(900))
	labels = [_mp_cell_label(c) for c in cells]
	top = min(top_labels, ns)
	yticklabels = nc == 0 ? string.(1:top) : [cells[1].labels[order[k]] for k in 1:top]
	ax2 = Axis(g[1, 1]; title = "(b) Sectoral total-import changes, % (" * clip_note_t * ")",
		xticks = nc == 0 ? ([1], [""]) : (1:nc, labels),
		yticks = (1:top, yticklabels),
		yreversed = true, xticklabelrotation = π / 4,
		xticklabelsize = 10, yticklabelsize = 9, titlesize = 14)
	hm = heatmap!(ax2, 1:max(nc, 1), 1:ns,
		nc == 0 ? fill(NaN, 1, ns) : mat;
		colormap = :RdBu, colorrange = cr)
	Colorbar(g[1, 2], hm; label = "%")
	exp_rel = [r for c in cells for r in c.exports if isfinite(r)]
	vdev = isempty(exp_rel) ? NaN : maximum(abs.(exp_rel .- 1.0))
	qdev = NaN
	if !isempty(cells)
		nb = length(ds.baseline.exports)
		trades = [(j <= nb && ds.baseline.exports[j] > 0.0) for j in 1:length(cells[1].exports)]
		q0 = let c0 = cells[1]
			Float64[(trades[j] && isfinite(e) && isfinite(pv) && pv != 0.0) ? e / pv : NaN
				for (j, (e, pv)) in enumerate(zip(c0.exports, c0.prices))]
		end
		best = 0.0
		found = false
		for c in cells
			length(c.exports) == length(c.prices) || continue
			length(c.exports) == length(q0) || continue
			for j in eachindex(q0)
				isfinite(q0[j]) && q0[j] != 0.0 || continue
				e = c.exports[j]
				pv = c.prices[j]
				(isfinite(e) && isfinite(pv) && pv != 0.0) || continue
				ratio = (e / pv) / q0[j]
				isfinite(ratio) || continue
				best = found ? max(best, abs(ratio - 1.0)) : abs(ratio - 1.0)
				found = true
			end
		end
		qdev = found ? best : NaN
	end
	trade_note = if !isfinite(vdev)
		"(c) Exports: no finite export relatives in this selection."
	else
		"(c) Exports are exogenous quantities; values move only with equilibrium prices " *
			"(max value deviation $(@sprintf("%.2f%%", 100.0 * vdev))" *
			(isfinite(qdev) ? "; max quantity deviation $(@sprintf("%.2e", qdev)) on trading sectors" : "") *
			")."
	end
	Label(fig[4, 1:2], trade_note, fontsize = 13)
	Label(fig[5, 1:2],
		_MP_BASELINE_NOTE * "  X = 100*gdp_x; M = -100*(gdp_m_final + gdp_m_int)." *
			"  Heatmap " * clip_note_t * ".",
		fontsize = 13)
	return fig
end

function BeyondHulten.save_matrix_figures(ds::BeyondHulten.MatrixDataset;
		outdir = "plots", prefix = ds.design, formats = ("png",), dpi = 150)
	mkpath(outdir)
	figs = [
		("overview", BeyondHulten.plot_matrix_overview(ds)),
		("wages", BeyondHulten.plot_matrix_wages(ds)),
		("prices", BeyondHulten.plot_matrix_prices(ds)),
		("quantities", BeyondHulten.plot_matrix_quantities(ds)),
		("consumption", BeyondHulten.plot_matrix_consumption(ds)),
		("trade", BeyondHulten.plot_matrix_trade(ds)),
	]
	ppu = clamp(round(Int, Float64(dpi) / 72.0), 1, 4)
	paths = String[]
	for (name, f) in figs, ext in formats
		p = joinpath(outdir, "$(prefix)_$(name).$(ext)")
		try
			save(p, f; px_per_unit = ppu)
		catch
			save(p, f)
		end
		push!(paths, p)
	end
	return paths
end
