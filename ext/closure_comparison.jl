# Physical units are points: the save helper maps 72 logical units to an inch.
function _cc_wrap(label, columns=24)
  lines, line = String[], ""
  for word in split(replace(label, r"[._]+" => " "))
    if !isempty(line) && length(line) + length(word) + 1 > columns
      push!(lines, line)
      line = ""
    end
    # Also wrap unbroken identifiers, preserving every character.
    while length(word) > columns
      !isempty(line) && push!(lines, line)
      line = ""
      push!(lines, first(word, columns))
      word = last(word, length(word) - columns)
    end
    line = isempty(line) ? string(word) : line * " " * word
  end
  !isempty(line) && push!(lines, line)
  return join(lines, '\n')
end

# Labour groups and financing styles follow the order of first appearance.
function _cc_group_labels(ids, labels, name)
  groups = unique(ids)
  labels === nothing && return groups, groups
  labels = string.(collect(labels))
  length(labels) == length(groups) || throw(DimensionMismatch(
    "$(name)_labels must have one entry per selected $name id"))
  all(s -> !isempty(strip(s)), labels) || throw(ArgumentError("display labels cannot be empty"))
  return groups, labels
end

const _CC_LINESTYLES = (:solid, :dash, :dot, :dashdot, :dashdotdot)

function OpenCGE.plot_closure_comparison(
  ds::OpenCGE.MatrixDataset,
  sectors::AbstractVector, closures::AbstractVector;
  financing=OpenCGE.matrix_financing_order(), sector_labels=nothing,
  closure_labels=nothing, labour_labels=nothing, financing_labels=nothing,
  width_cm::Real=17, axis_height_cm::Real=7, fontsize::Real=10
)
  isfinite(width_cm) && width_cm >= 12 || throw(ArgumentError("width_cm must be finite and at least 12 cm"))
  isfinite(axis_height_cm) && axis_height_cm >= 3 || throw(ArgumentError("axis_height_cm must be finite and at least 3 cm"))
  isfinite(fontsize) && fontsize >= 8 || throw(ArgumentError("fontsize must be finite and at least 8 print points"))
  data = OpenCGE.closure_comparison_data(
    ds, sectors, closures;
    financing, sector_labels, closure_labels
  )
  labours, llabels = _cc_group_labels(data.labours, labour_labels, "labour")
  fins, flabels = _cc_group_labels(data.financings, financing_labels, "financing")
  length(fins) <= length(_CC_LINESTYLES) || throw(ArgumentError(
    "at most $(length(_CC_LINESTYLES)) financing ids can be distinguished; got $(length(fins))"))
  ns, nl = size(data.changes, 1), length(labours)
  members = [findall(==(l), data.labours) for l in labours]
  fin_of = [findfirst(==(f), fins) for f in data.financings]

  # Bars sit on unit slots; groups and sectors are separated by wider gaps.
  slot, barwidth, groupgap, sectorgap = 1.0, 0.8, 0.6, 1.4
  x = zeros(ns, nl, 3)
  pos = 0.0
  for i in 1:ns, g in 1:nl, v in 1:3
    pos += slot + (v == 1 && g > 1 ? groupgap : 0.0) +
      (v == 1 && g == 1 && i > 1 ? sectorgap : 0.0)
    x[i, g, v] = pos
  end
  lower = [minimum(data.changes[i, v, k] for k in members[g]) for i in 1:ns, g in 1:nl, v in 1:3]
  upper = [maximum(data.changes[i, v, k] for k in members[g]) for i in 1:ns, g in 1:nl, v in 1:3]

  width = round(Int, width_cm / 2.54 * 72)
  axis_height = axis_height_cm / 2.54 * 72
  # Approximate plot width (minus y ticks and padding) sets label wrapping.
  pt_per_unit = (width - 70) / (pos + slot)
  group_columns = max(3, floor(Int, 3 * slot * pt_per_unit / (0.6 * fontsize)))
  sector_columns = max(3, floor(Int, (3nl * slot + (nl - 1) * groupgap) * pt_per_unit / (0.6 * fontsize)))
  fig = Figure(
    size=(width, round(Int, axis_height + 160)), fontsize=fontsize,
    figure_padding=(10, 10, 8, 8), backgroundcolor=:white
  )
  colors = Makie.wong_colors()[[2, 3, 4]]
  Legend(
    fig[1, 1], [PolyElement(color=(c, 0.28), strokewidth=0) for c in colors],
    [L"\Delta y", L"\Delta W_D", L"\Delta p"];
    orientation=:horizontal, framevisible=false, labelsize=fontsize,
    patchsize=(14, 10), colgap=12, padding=(0, 0, 0, 0), tellwidth=false
  )
  financing_legend = Legend(
    fig[2, 1], [
      LineElement(color=:gray20, linestyle=_CC_LINESTYLES[f], linewidth=1.3)
      for f in eachindex(fins)
    ], _cc_wrap.(flabels, 26);
    orientation=:horizontal, framevisible=false, labelsize=fontsize,
    patchsize=(20, 10), colgap=12, padding=(0, 0, 0, 0), tellwidth=false
  )
  # Use the measured legend width for unusually wide labels or large text.
  if financing_legend.layoutobservables.computedbbox[].widths[1] > width - 20
    financing_legend.orientation[] = :vertical
    financing_legend.nbanks[] = 2
  end
  group_centres = vec(permutedims(dropdims(sum(x; dims=3); dims=3) ./ 3))
  sector_centres = [(x[i, 1, 1] + x[i, nl, 3]) / 2 for i in 1:ns]
  ax = Axis(
    fig[3, 1]; ylabel="Change vs. baseline (%)", height=axis_height,
    xticks=(group_centres, repeat(_cc_wrap.(llabels, group_columns), ns)),
    xgridvisible=false, xticksvisible=false,
    xticklabelsize=fontsize, yticklabelsize=fontsize,
    xlabelsize=fontsize, ylabelsize=fontsize,
    ygridcolor=(:black, 0.08), topspinevisible=false, rightspinevisible=false
  )
  # Sector names sit above the plot on a linked, otherwise empty axis.
  sector_axis = Axis(
    fig[3, 1]; height=axis_height, xaxisposition=:top,
    xticks=(sector_centres, _cc_wrap.(data.labels, sector_columns)),
    xticksvisible=false, xticklabelsize=fontsize, backgroundcolor=:transparent,
    xgridvisible=false, ygridvisible=false, yticklabelsvisible=false,
    yticksvisible=false, leftspinevisible=false, rightspinevisible=false,
    topspinevisible=false, bottomspinevisible=false
  )
  linkaxes!(ax, sector_axis)
  hlines!(ax, [0.0]; color=:gray50, linewidth=0.8, linestyle=:dash)
  ns > 1 && vlines!(
    ax, [(x[i, nl, 3] + x[i+1, 1, 1]) / 2 for i in 1:(ns-1)];
    color=(:black, 0.1), linewidth=0.75
  )
  for v in 1:3
    barplot!(
      ax, vec(x[:, :, v]'), vec(upper[:, :, v]'); fillto=vec(lower[:, :, v]'),
      width=barwidth, gap=0, color=(colors[v], 0.28), strokewidth=0
    )
  end
  # One line per financing cell; lines also mark zero-height ranges.
  for f in eachindex(fins)
    segments, segment_colors = Point2f[], RGBAf[]
    for i in 1:ns, g in 1:nl, k in members[g], v in 1:3
      fin_of[k] == f || continue
      y = data.changes[i, v, k]
      push!(
        segments, Point2f(x[i, g, v] - barwidth / 2, y),
        Point2f(x[i, g, v] + barwidth / 2, y)
      )
      append!(segment_colors, fill(Makie.to_color(colors[v]), 2))
    end
    linesegments!(
      ax, segments; color=segment_colors,
      linestyle=_CC_LINESTYLES[f], linewidth=1.3
    )
  end
  values = vcat(0.0, vec(lower), vec(upper))
  lo, hi, _ = _mp_panel_span(
    values;
    min_span=all(iszero, values) ? 1.0 : 0.0, pad=0.08
  )
  xlims!(ax, 0.5 * slot, pos + 0.5 * slot)
  ylims!(ax, lo, hi)
  ax.yticks = Makie.WilkinsonTicks(5; k_min=4, k_max=6)
  # Match the actual tick locations, including half-percent steps at wide
  # spans, and preserve tiny effects with scientific notation.
  ticks = Makie.get_tickvalues(ax.yticks[], lo, hi)
  digits = 0
  while digits < 8 && any(t -> abs(t - round(t; digits)) > (hi - lo) * 1.0e-8, ticks)
    digits += 1
  end
  ax.ytickformat = hi - lo < 1.0e-4 || hi - lo > 1.0e5 ? "{:.2e}%" : "{:.$(digits)f}%"
  rowgap!(fig.layout, 6)
  Makie.resize_to_layout!(fig)
  return fig
end

function OpenCGE.save_closure_comparison(
  path::AbstractString,
  ds::OpenCGE.MatrixDataset, sectors::AbstractVector, closures::AbstractVector;
  dpi::Real=600, kwargs...
)
  isfinite(dpi) && dpi >= 300 || throw(ArgumentError("dpi must be finite and at least 300 for print"))
  lowercase(splitext(path)[2]) == ".png" || throw(ArgumentError("save_closure_comparison writes .png files"))
  fig = OpenCGE.plot_closure_comparison(ds, sectors, closures; kwargs...)
  # FileIO forwards dpi to PNGFiles, which writes the PNG pHYs metadata.
  img = Makie.colorbuffer(fig; px_per_unit=dpi / 72)
  GLMakie.FileIO.save(path, img; dpi)
  return path
end
