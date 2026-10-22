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
    clabels = _cc_wrap.(data.closure_labels, 26)
    width = round(Int, width_cm / 2.54 * 72)
    # A separate lane for every closure keeps coincident responses readable.
    barheight = max(9.0, nc * 3.0)
    rowheight = max(54.0, 3 * (barheight + 4) + 10,
        maximum(count(==('\n'), l) + 1 for l in labels) * fontsize * 1.2 + 6)
    legend_columns = max(1, min(3, floor(Int, width / (fontsize * 15))))
    legend_rows = cld(nc, legend_columns)
    keyheight = max(12, nc * 3)
    legend_height = legend_rows * (max(fontsize * 1.6, keyheight) + 3) + fontsize * 4.8 +
        sum(count(==('\n'), l) for l in clabels) * fontsize * 1.2
    height = ceil(Int, ns * rowheight + legend_height + fontsize * 5)
    fig = Figure(size = (width, height), fontsize = fontsize,
        figure_padding = (10, 10, 8, 8), backgroundcolor = :white)
    colors = Makie.wong_colors()[[2, 3, 4]]
    Legend(fig[1, 1], [PolyElement(polycolor = c) for c in colors],
        ["Output", "Consumption (W_D)", "Price"];
        orientation = :horizontal, framevisible = false, labelsize = fontsize,
        patchsize = (14, 10), colgap = 12, padding = (0, 0, 0, 0), tellwidth = false)
    # Lane location in each key matches the lane inside every coloured bar.
    elements = [LineElement(color = :black, linewidth = 0.8,
        linepoints = Point2f[(0.5, 1 - (k - 0.85) / nc),
            (0.5, 1 - (k - 0.15) / nc)]) for k in 1:nc]
    Legend(fig[2, 1], elements, clabels, "Closure ticks (top to bottom)";
        orientation = :vertical, nbanks = legend_columns, tellheight = true,
        framevisible = false, labelsize = fontsize, titlesize = fontsize,
        patchsize = (12, keyheight), colgap = 12, rowgap = 3,
        padding = (0, 0, 0, 0), tellwidth = false)
    ax = Axis(fig[3, 1]; xlabel = "Change vs. baseline (%)",
        yticks = (1:ns, labels), yreversed = true, ygridvisible = false,
        xticklabelsize = fontsize, yticklabelsize = fontsize,
        xlabelsize = fontsize, ylabelsize = fontsize,
        xgridcolor = (:black, 0.08), topspinevisible = false, rightspinevisible = false)
    vlines!(ax, [0.0]; color = :gray50, linewidth = 0.8, linestyle = :dash)
    band = barheight / rowheight
    spacing = (barheight + 4) / rowheight
    for v in 1:3
        positions = collect(1.0:ns) .+ (v - 2) * spacing
        barplot!(ax, positions, data.maxima[:, v]; fillto = data.minima[:, v],
            direction = :x, width = band, gap = 0, color = colors[v],
            strokecolor = colors[v], strokewidth = 0.5)
        for k in 1:nc
            points = Point2f[]
            for i in 1:ns
                center = positions[i] - band / 2 + (k - 0.5) * band / nc
                half = 0.35 * band / nc
                x = data.changes[i, v, k]
                push!(points, Point2f(x, center - half), Point2f(x, center + half))
            end
            # Keep ticks in front of the filled bars, including at tied extrema.
            linesegments!(ax, points; color = :black, linewidth = 0.8,
                depth_shift = -1f-3)
        end
    end
    values = vcat(0.0, vec(data.minima), vec(data.maxima))
    lo, hi, _ = _mp_panel_span(values;
        min_span = all(iszero, values) ? 1.0 : 0.0, pad = 0.08)
    xlims!(ax, lo, hi)
    ylims!(ax, 0.5, ns + 0.5)
    ax.yreversed = true
    ax.xticks = Makie.WilkinsonTicks(5; k_min = 4, k_max = 6)
    # Preserve tiny effects instead of rounding every tick to zero.
    digits = clamp(ceil(Int, -log10((hi - lo) / 5)), 0, 8)
    ax.xtickformat = hi - lo < 1e-4 || hi - lo > 1e5 ? "{:.2e}%" : "{:.$(digits)f}%"
    rowgap!(fig.layout, 6)
    Makie.update_state_before_display!(fig)
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
