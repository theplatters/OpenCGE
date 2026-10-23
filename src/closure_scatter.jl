"""
    closure_scatter_data(ds, closure_pairs; programme, domar_weights = nothing,
        sector_labels = nothing, closure_labels = nothing)

Prepare output and price percentage changes for [`plot_closure_scatter`](@ref).
`closure_pairs` selects explicit `(labour, financing)` tuples or `labour => financing`
Pairs, using strings or symbols, in the requested order. Each pair must match
exactly one `MatrixDataset` cell; missing, duplicate, or ambiguous pairs are errors.

Attach a top-legend name with `"Benchmark" => (:BF, :F1)`. Named and unnamed
entries may be mixed; unnamed entries use the labour–financing ids. Optional
`closure_labels` overrides all attached names in selection order. Names must
be nonempty and unique across the selection.

All sectors are retained in dataset order. Every nonzero `programme` entry
identifies a directly shocked sector, including negative shocks. Returns
`output[sector, closure]`, `price[sector, closure]`, `shocked`, `domar_weights`,
`closure_pairs`, and identifying labels. Nonfinite plotted responses are errors.

`domar_weights` optionally supplies pre-shock sectoral Domar weights (e.g.
`data.λ`). By default, baseline `prices .* quantities` supplies gross-output
values proportional to those weights: the common baseline GDP denominator
cancels when sizing points. Weights must be finite and nonnegative, with at
least one positive entry. They are shared by all closures and shock groups.
"""
function closure_scatter_data(ds::MatrixDataset, closure_pairs::AbstractVector;
        programme::AbstractVector, domar_weights = nothing,
        sector_labels = nothing, closure_labels = nothing)
    isempty(ds.cells) && throw(ArgumentError("dataset contains no cells"))
    isempty(closure_pairs) && throw(ArgumentError("select at least one labour–financing pair"))
    pairs = Tuple{String,String}[]
    pair_labels = String[]
    run_ids = String[]
    for entry in closure_pairs
        pair, label = entry, nothing
        if entry isa Pair && (last(entry) isa Tuple || last(entry) isa Pair)
            label, pair = first(entry), last(entry)
            label isa AbstractString || label isa Symbol || throw(ArgumentError(
                "an attached closure label must be a string or symbol"))
        end
        (pair isa Pair || (pair isa Tuple && length(pair) == 2)) || throw(ArgumentError(
            "closures must be explicit (labour, financing) tuples or labour => financing pairs"))
        labour, financing = first(pair), last(pair)
        all(x -> x isa AbstractString || x isa Symbol, (labour, financing)) ||
            throw(ArgumentError("labour and financing ids must be strings or symbols"))
        key = (string(labour), string(financing))
        key in pairs && throw(ArgumentError("closure pair $key was selected twice"))
        matches = filter(c -> (c.labour, c.financing) == key, ds.cells)
        length(matches) == 1 || throw(ArgumentError(
            "expected one cell for $(key[1])/$(key[2]), found $(length(matches)); select a dataset with one variant per pair"))
        push!(pairs, key)
        push!(pair_labels, label === nothing ? key[1] * "–" * key[2] : string(label))
        push!(run_ids, only(matches).run_id)
    end
    n = length(first(ds.cells).sectors)
    (; ref, cells, labels, clabels) = _closure_comparison_selection(
        ds, collect(1:n), run_ids; financing = String[], sector_labels,
        closure_labels = closure_labels === nothing ? pair_labels : closure_labels)
    length(programme) == n || throw(DimensionMismatch("programme must have $n entries"))
    all(isfinite, programme) || throw(ArgumentError("programme must contain only finite entries"))
    weights = if domar_weights === nothing
        length(ds.baseline.prices) == length(ds.baseline.quantities) == n ||
            throw(DimensionMismatch("baseline prices and quantities must have $n entries"))
        all(x -> isfinite(x) && x >= 0, ds.baseline.prices) &&
            all(x -> isfinite(x) && x >= 0, ds.baseline.quantities) ||
            throw(ArgumentError("baseline prices and quantities must be finite and nonnegative"))
        ds.baseline.prices .* ds.baseline.quantities
    else
        Float64.(collect(domar_weights))
    end
    length(weights) == n || throw(DimensionMismatch("domar_weights must have $n entries"))
    all(w -> isfinite(w) && w >= 0, weights) && any(>(0), weights) ||
        throw(ArgumentError("Domar weights must be finite and nonnegative, with at least one positive entry"))
    output, price = zeros(n, length(cells)), zeros(n, length(cells))
    for (k, cell) in enumerate(cells)
        cell.sectors == ref.sectors || throw(ArgumentError("sector order differs in $(cell.run_id)"))
        for (values, dest, variable) in ((cell.quantities, output, "output"),
                (cell.prices, price, "price"))
            length(values) == n || throw(DimensionMismatch("$variable vector in $(cell.run_id) must have $n entries"))
            for j in 1:n
                change = 100.0 * (values[j] - 1.0)
                isfinite(change) || throw(ArgumentError(
                    "undefined/nonfinite $variable in $(cell.run_id), sector $(ref.sectors[j])"))
                dest[j, k] = change
            end
        end
    end
    return (; sectors = copy(ref.sectors), labels, run_ids, closure_labels = clabels,
        closure_pairs = pairs, domar_weights = weights, output, price,
        shocked = .!iszero.(programme))
end

"""
    plot_closure_scatter(ds, closure_pairs; programme, domar_weights = nothing,
        sector_labels = nothing, closure_labels = nothing, closure_colors = nothing,
        unshocked_alpha = 0.35, max_markersize = 24,
        size = (1100, 650), fontsize = 16)

Compare explicit labour–financing closure pairs on one output-change (%) ×
price-change (%) scatter. Every sector is a circular point; each closure pair
has its own color. Point **area** is proportional to the pre-shock Domar weight,
with the largest sector having diameter `max_markersize` in pixels. The same
sector has the same size under every closure, regardless of shock status.
Shocked sectors use darker opaque shades of the closure color; other sectors
use lighter shades at `unshocked_alpha` opacity and are drawn first.
Coordinates are never jittered.

Inputs follow [`closure_scatter_data`](@ref). Default colors are stable when
reordering or subsetting pairs in the same dataset. `closure_colors` optionally
supplies one base color per selected pair. Attach top-legend names directly:
`["Benchmark" => (:BF, :F1), "Alternative" => (:BETA, :F2)]`.
The `closure_labels` keyword can also supply or override all names.
Requires `using GLMakie`; returns a
Makie `Figure`, which can be saved with `save("closures.png", fig)`.
"""
function plot_closure_scatter(args...; kwargs...)
    _require_glmakie(:plot_closure_scatter)
end
