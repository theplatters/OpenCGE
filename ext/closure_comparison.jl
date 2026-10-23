# Physical units are points: the save helper maps 72 logical units to an inch.
function _cc_wrap(label, columns = 24)
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

# Only nearby responses move off the centreline. Keep even large tie groups
# inside the fixed band, preserving the selected closure order within a tie.
function _cc_tie_offsets(values, tolerance, barheight)
    offsets = zeros(length(values))
    order = sortperm(values)
    start = 1
    while start <= length(order)
        stop = start
        while stop < length(order) && values[order[stop + 1]] - values[order[stop]] <= tolerance
            stop += 1
        end
        n = stop - start + 1
        if n > 1
            step = min(2.2, (barheight - 1) / (n - 1))
            for (j, k) in enumerate(sort(order[start:stop]))
                offsets[k] = (j - (n + 1) / 2) * step
            end
        end
        start = stop + 1
    end
    return offsets
end

function OpenCGE.plot_closure_comparison(ds::OpenCGE.MatrixDataset,
        sectors::AbstractVector, closures::AbstractVector;
        financing = OpenCGE.matrix_financing_order(), sector_labels = nothing,
        closure_labels = nothing, width_cm::Real = 17, fontsize::Real = 10)
    isfinite(width_cm) && width_cm >= 12 || throw(ArgumentError("width_cm must be finite and at least 12 cm"))
    isfinite(fontsize) && fontsize >= 8 || throw(ArgumentError("fontsize must be finite and at least 8 print points"))
    data = OpenCGE.closure_comparison_data(ds, sectors, closures;
        financing, sector_labels, closure_labels)
    ns, _, nc = size(data.changes)
    labels = _cc_wrap.(data.labels)
    width = round(Int, width_cm / 2.54 * 72)
    barheight = 9.0
    rowheight = max(56.0,
        maximum(count(==('\n'), l) + 1 for l in labels) * fontsize * 1.2 + 6)
    legend_columns = min(nc, maximum(length, data.closure_labels) > 16 ? 3 : 4)
    wrap_columns = clamp(floor(Int,
        ((width - 20 - 12 * (legend_columns - 1)) / legend_columns - 20) / (0.6 * fontsize)), 1, 26)
    clabels = _cc_wrap.(data.closure_labels, wrap_columns)
    height = ceil(Int, ns * rowheight + 100)
    fig = Figure(size = (width, height), fontsize = fontsize,
        figure_padding = (10, 10, 8, 8), backgroundcolor = :white)
    colors = Makie.wong_colors()[[2, 3, 4]]
    Legend(fig[1, 1], [PolyElement(color = (c, 0.28),
            strokecolor = (c, 0.85), strokewidth = 1.0) for c in colors],
        ["Output", "Consumption (W_D)", "Price"];
        orientation = :horizontal, framevisible = false, labelsize = fontsize,
        patchsize = (14, 10), colgap = 12, padding = (0, 0, 0, 0), tellwidth = false)
    shapes = [:circle, :rect, :utriangle, :diamond, :dtriangle,
        :pentagon, :xcross, :star5, :hexagon]
    markers = [shapes[mod1(k, length(shapes))] for k in 1:nc]
    fills = [isodd(cld(k, length(shapes))) ? :gray20 : :white for k in 1:nc]
    elements = [MarkerElement(marker = markers[k], color = fills[k],
        strokecolor = :gray20, markersize = 7.5, strokewidth = 0.9) for k in 1:nc]
    closure_legend = Legend(fig[2, 1], elements, clabels;
        orientation = :vertical, nbanks = legend_columns, tellheight = true,
        framevisible = false, labelsize = fontsize, titlesize = fontsize,
        patchsize = (14, 12), colgap = 12, rowgap = 2,
        padding = (0, 0, 0, 0), tellwidth = false)
    # Use the measured legend width for unusually wide glyphs or large text.
    while closure_legend.layoutobservables.computedbbox[].widths[1] > width - 20 &&
            closure_legend.nbanks[] > 1
        closure_legend.nbanks[] -= 1
    end
    ax = Axis(fig[3, 1]; xlabel = "Change vs. baseline (%)",
        height = ns * rowheight,
        yticks = (1:ns, labels), yreversed = true, ygridvisible = false,
        xticklabelsize = fontsize, yticklabelsize = fontsize,
        xlabelsize = fontsize, ylabelsize = fontsize,
        xgridcolor = (:black, 0.08), topspinevisible = false, rightspinevisible = false)
    vlines!(ax, [0.0]; color = :gray50, linewidth = 0.8, linestyle = :dash)
    hlines!(ax, collect(1.5:(ns - 0.5)); color = (:black, 0.10), linewidth = 0.75)
    values = vcat(0.0, vec(data.minima), vec(data.maxima))
    lo, hi, _ = _mp_panel_span(values;
        min_span = all(iszero, values) ? 1.0 : 0.0, pad = 0.08)
    band = barheight / rowheight
    spacing = 14.0 / rowheight
    for v in 1:3
        positions = collect(1.0:ns) .+ (v - 2) * spacing
        barplot!(ax, positions, data.maxima[:, v]; fillto = data.minima[:, v],
            direction = :x, width = band, gap = 0, color = (colors[v], 0.28),
            strokecolor = (colors[v], 0.85), strokewidth = 1.0)
        caps = Point2f[]
        offsets = zeros(ns, nc)
        for i in 1:ns
            for x in unique((data.minima[i, v], data.maxima[i, v]))
                push!(caps, Point2f(x, positions[i] - band / 2),
                    Point2f(x, positions[i] + band / 2))
            end
            offsets[i, :] = _cc_tie_offsets(view(data.changes, i, v, :),
                0.015 * (hi - lo), barheight) ./ rowheight
        end
        # Explicit caps also keep zero-width ranges visible.
        linesegments!(ax, caps; color = colors[v], linewidth = 1.0)
        for k in 1:nc
            scatter!(ax, data.changes[:, v, k], positions .+ offsets[:, k];
                marker = markers[k], color = fills[k], strokecolor = :gray20,
                markersize = 7.5, strokewidth = 0.9,
                depth_shift = -1f-3)
        end
    end
    xlims!(ax, lo, hi)
    ylims!(ax, 0.5, ns + 0.5)
    ax.yreversed = true
    ax.xticks = Makie.WilkinsonTicks(5; k_min = 4, k_max = 6)
    # Match the actual tick locations, including half-percent steps at wide
    # spans, and preserve tiny effects with scientific notation.
    ticks = Makie.get_tickvalues(ax.xticks[], lo, hi)
    digits = 0
    while digits < 8 && any(t -> abs(t - round(t; digits)) > (hi - lo) * 1e-8, ticks)
        digits += 1
    end
    ax.xtickformat = hi - lo < 1e-4 || hi - lo > 1e5 ? "{:.2e}%" : "{:.$(digits)f}%"
    rowgap!(fig.layout, 6)
    Makie.resize_to_layout!(fig)
    return fig
end

function OpenCGE.save_closure_comparison(path::AbstractString,
        ds::OpenCGE.MatrixDataset, sectors::AbstractVector, closures::AbstractVector;
        dpi::Real = 600, kwargs...)
    isfinite(dpi) && dpi >= 300 || throw(ArgumentError("dpi must be finite and at least 300 for print"))
    lowercase(splitext(path)[2]) == ".png" || throw(ArgumentError("save_closure_comparison writes .png files"))
    fig = OpenCGE.plot_closure_comparison(ds, sectors, closures; kwargs...)
    # FileIO forwards dpi to PNGFiles, which writes the PNG pHYs metadata.
    img = Makie.colorbuffer(fig; px_per_unit = dpi / 72)
    GLMakie.FileIO.save(path, img; dpi)
    return path
end
