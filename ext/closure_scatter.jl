# Derive precision from actual ticks, so 2.5% never appears as 3%.
function _cs_tickformat(ticks, lo, hi)
    if 1e-4 <= hi - lo <= 1e5
        values = Makie.get_tickvalues(ticks, lo, hi)
        digits = 0
        while digits < 8 && any(t -> abs(t - round(t; digits)) > (hi - lo) * 1e-8, values)
            digits += 1
        end
        return "{:.$(digits)f}%"
    end
    return "{:.2e}%"
end

function OpenCGE.plot_closure_scatter(ds::OpenCGE.MatrixDataset,
        closure_pairs::AbstractVector; programme::AbstractVector, domar_weights = nothing,
        sector_labels = nothing, closure_labels = nothing, closure_colors = nothing,
        unshocked_alpha::Real = 0.35, max_markersize::Real = 24,
        size = (1100, 650), fontsize::Real = 16)
    isfinite(unshocked_alpha) && 0 < unshocked_alpha < 1 || throw(ArgumentError(
        "unshocked_alpha must be between 0 and 1 (exclusive)"))
    isfinite(max_markersize) && max_markersize > 0 || throw(ArgumentError(
        "max_markersize must be finite and positive"))
    isfinite(fontsize) && fontsize > 0 || throw(ArgumentError("fontsize must be finite and positive"))
    data = OpenCGE.closure_scatter_data(ds, closure_pairs;
        programme, domar_weights, sector_labels, closure_labels)
    _, nc = Base.size(data.output)
    base = if closure_colors === nothing
        # Assign colors using the full dataset so pair colors survive reordering
        # and subsetting. The selected pairs are never expanded to a cross product.
        all_pairs = sort!(unique([(c.labour, c.financing) for c in ds.cells]))
        palette = length(all_pairs) <= 7 ? Makie.wong_colors()[1:length(all_pairs)] :
            Makie.distinguishable_colors(length(all_pairs),
                [RGBf(1, 1, 1), RGBf(0, 0, 0)]; dropseed = true)
        [palette[findfirst(==(pair), all_pairs)] for pair in data.closure_pairs]
    else
        closure_colors
    end
    length(base) == nc || throw(DimensionMismatch("closure_colors must have one entry per closure pair"))
    rgb = RGBf.(Makie.to_color.(base))
    length(unique(rgb)) == nc || throw(ArgumentError("closure_colors must be distinct"))
    dark = [RGBAf(0.75c.r, 0.75c.g, 0.75c.b, 1) for c in rgb]
    light = [RGBAf(0.65c.r + 0.35, 0.65c.g + 0.35, 0.65c.b + 0.35,
        unshocked_alpha) for c in rgb]
    # Makie markersize is diameter: square roots preserve relative point areas.
    diameters = max_markersize .* sqrt.(data.domar_weights ./ maximum(data.domar_weights))
    fig = Figure(; size, fontsize, backgroundcolor = :white)
    Legend(fig[1, 1], [MarkerElement(marker = :circle, color = dark[k],
            markersize = 12) for k in 1:nc], _cc_wrap.(data.closure_labels), "Labour–financing closures";
        orientation = :vertical, nbanks = min(nc, 4), framevisible = false,
        tellwidth = false, tellheight = true, labelsize = fontsize, titlesize = fontsize)
    ax = Axis(fig[2, 1]; xlabel = "Output change vs. baseline (%)",
        ylabel = "Price change vs. baseline (%)", xlabelsize = fontsize,
        ylabelsize = fontsize, xticklabelsize = fontsize, yticklabelsize = fontsize,
        topspinevisible = false, rightspinevisible = false)
    hlines!(ax, [0.0]; color = :gray60, linestyle = :dash, linewidth = 1)
    vlines!(ax, [0.0]; color = :gray60, linestyle = :dash, linewidth = 1)
    # Draw every background closure before shocked sectors; only shade and alpha
    # encode shock status, so point areas always retain the Domar-weight ratio.
    for shocked in (false, true)
        idx = findall(==(shocked), data.shocked)
        isempty(idx) && continue
        for k in 1:nc
            scatter!(ax, data.output[idx, k], data.price[idx, k];
                marker = :circle, color = shocked ? dark[k] : light[k],
                markersize = diameters[idx], strokewidth = 0)
        end
    end
    xlo, xhi, _ = _mp_panel_span(vcat(0.0, vec(data.output)); pad = 0.10)
    ylo, yhi, _ = _mp_panel_span(vcat(0.0, vec(data.price)); pad = 0.10)
    xlims!(ax, xlo, xhi)
    ylims!(ax, ylo, yhi)
    ticks = Makie.WilkinsonTicks(5; k_min = 4, k_max = 6)
    ax.xticks = ticks
    ax.yticks = ticks
    ax.xtickformat = _cs_tickformat(ticks, xlo, xhi)
    ax.ytickformat = _cs_tickformat(ticks, ylo, yhi)
    Legend(fig[3, 1],
        [MarkerElement(marker = :circle, color = :gray25, markersize = 12),
            MarkerElement(marker = :circle, color = (:gray65, unshocked_alpha), markersize = 12)],
        ["Shocked sector", "Other sector"], "Point area ∝ pre-shock Domar weight";
        orientation = :horizontal, framevisible = false, labelsize = fontsize,
        titlesize = fontsize - 1, tellwidth = false, tellheight = true)
    rowgap!(fig.layout, 12)
    return fig
end
