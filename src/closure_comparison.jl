"""
    closure_comparison_data(ds, sectors, closures;
        financing = matrix_financing_order(), sector_labels = nothing,
        closure_labels = nothing)

Prepare sectoral percentage changes for [`plot_closure_comparison`](@ref),
without loading a plotting backend. `ds` is a `MatrixDataset`, whose cell
vectors already contain level / baseline ratios. Consumption means sectoral
household consumption, not aggregate welfare.

`sectors` selects integer indices, exact sector identifiers, or exact display
labels, in the requested order. `sector_labels` optionally supplies one curated
label per selected sector. `closures` selects exact run ids or labour ids
(strings or symbols). A labour id expands to every requested `financing` id,
in financing order; missing or ambiguous cells are errors. Use exact run ids
to select a subset or distinguish parameter variants. `closure_labels` supplies
one label per expanded cell. No solutions are computed or silently omitted.

Returns `changes[sector, variable, closure]`, `minima[sector, variable]`,
`maxima[sector, variable]`, identifying labels, and the `labours` and
`financings` ids of each expanded cell. Variables are output,
household consumption, and price, in that order. Ranges are the exact extrema
over the selected cells; zero is not included unless a cell attains it.
Undefined/nonfinite selected observations are rejected.
"""
function closure_comparison_data(ds::MatrixDataset, sectors::AbstractVector,
        closures::AbstractVector; financing = matrix_financing_order(),
        sector_labels = nothing, closure_labels = nothing)
    (; ref, idx, cells, labels, clabels, run_ids) = _closure_comparison_selection(
        ds, sectors, closures; financing, sector_labels, closure_labels)
    n = length(ref.sectors)
    changes = Array{Float64}(undef, length(idx), 3, length(cells))
    for (k, cell) in enumerate(cells)
        cell.sectors == ref.sectors || throw(ArgumentError("sector order differs in $(cell.run_id)"))
        for (v, values) in enumerate((cell.quantities, cell.consumption, cell.prices))
            length(values) == n || throw(DimensionMismatch("sector vector in $(cell.run_id) must have $n entries"))
            for (i, j) in enumerate(idx)
                change = 100.0 * (values[j] - 1.0)
                isfinite(change) || throw(ArgumentError(
                    "undefined/nonfinite response in $(cell.run_id), sector $(ref.sectors[j]), variable $v"))
                changes[i, v, k] = change
            end
        end
    end
    return (; sectors = ref.sectors[idx], labels, run_ids, closure_labels = clabels,
        labours = [c.labour for c in cells], financings = [c.financing for c in cells],
        variables = (:output, :consumption, :price), changes,
        minima = dropdims(minimum(changes; dims = 3); dims = 3),
        maxima = dropdims(maximum(changes; dims = 3); dims = 3))
end

# Shared selectors keep the range and scatter plots on the same cell ordering.
function _closure_comparison_selection(ds::MatrixDataset, sectors::AbstractVector,
        closures::AbstractVector; financing, sector_labels, closure_labels)
    isempty(sectors) && throw(ArgumentError("select at least one sector"))
    isempty(closures) && throw(ArgumentError("select at least one closure"))
    isempty(ds.cells) && throw(ArgumentError("dataset contains no cells"))
    ref = first(ds.cells)
    n = length(ref.sectors)
    length(ref.labels) == n || throw(DimensionMismatch("sector identifiers and labels must match"))
    length(unique(ref.sectors)) == n || throw(ArgumentError("sector identifiers must be unique"))
    idx = Int[]
    for sector in sectors
        j = if sector isa Integer && !(sector isa Bool)
            1 <= sector <= n || throw(ArgumentError("sector index $sector is outside 1:$n"))
            Int(sector)
        elseif sector isa AbstractString || sector isa Symbol
            key = string(sector)
            matches = findall(==(key), ref.sectors)
            isempty(matches) && (matches = findall(==(key), ref.labels))
            length(matches) == 1 || throw(ArgumentError(
                "sector $(repr(key)) is unknown or ambiguous; use an integer index"))
            only(matches)
        else
            throw(ArgumentError("sectors must be integer indices or sector names"))
        end
        j in idx && throw(ArgumentError("sector $(ref.sectors[j]) was selected twice"))
        push!(idx, j)
    end
    fins = string.(collect(financing))
    length(unique(fins)) == length(fins) || throw(ArgumentError("financing ids must be unique"))
    cells = MatrixCellData[]
    for closure in closures
        key = string(closure)
        matches = filter(c -> c.run_id == key, ds.cells)
        if isempty(matches)
            isempty(fins) && throw(ArgumentError("financing cannot be empty when selecting labour ids"))
            for fin in fins
                candidates = filter(c -> c.labour == key && c.financing == fin, ds.cells)
                length(candidates) == 1 || throw(ArgumentError(
                    "expected one cell for $key/$fin, found $(length(candidates)); use exact run ids for variants"))
                push!(cells, only(candidates))
            end
        else
            length(matches) == 1 || throw(ArgumentError("run id $key is not unique"))
            push!(cells, only(matches))
        end
    end
    run_ids = [c.run_id for c in cells]
    length(unique(run_ids)) == length(run_ids) || throw(ArgumentError("a closure cell was selected twice"))
    labels = sector_labels === nothing ? ref.labels[idx] : string.(collect(sector_labels))
    length(labels) == length(idx) || throw(DimensionMismatch("sector_labels must have one entry per selected sector"))
    clabels = closure_labels === nothing ? [c.labour * "–" * c.financing for c in cells] :
        string.(collect(closure_labels))
    length(clabels) == length(cells) || throw(DimensionMismatch("closure_labels must have one entry per expanded closure cell"))
    length(unique(clabels)) == length(clabels) || throw(ArgumentError(
        "closure labels must be unique; supply closure_labels to distinguish parameter variants"))
    all(s -> !isempty(strip(s)), vcat(labels, clabels)) || throw(ArgumentError("display labels cannot be empty"))
    return (; ref, idx, cells, labels, clabels, run_ids)
end

"""
    plot_closure_comparison(ds, sectors, closures; width_cm = 17,
        axis_height_cm = 7, fontsize = 10, financing = matrix_financing_order(),
        sector_labels = nothing, closure_labels = nothing,
        labour_labels = nothing, financing_labels = nothing)

Publication-sized vertical range bars. Each selected sector holds one group of
bars per labour closure, and each group holds output, household consumption
(W_D), and price bars in the existing panel palette. A translucent bar runs
from the smallest to the largest result over that labour closure's selected
financing cells; a line across the bar marks each financing cell's result,
with one line style per financing id (at most five). Labour groups and
financing styles follow the order of first appearance in the selection. All
three variables are retained, including unchanged prices. The dashed
horizontal reference is zero change.

Selectors follow [`closure_comparison_data`](@ref); `closure_labels` is
validated but not drawn. `labour_labels` and `financing_labels` supply one
display label per distinct labour or financing id. Font sizes are in print
points at `width_cm`; `axis_height_cm` is the plot-area height. Long labels
wrap without truncation. For long selections, split sectors across figures.
Requires `using GLMakie`. Returns a Makie `Figure`; see
[`save_closure_comparison`](@ref) for a PNG at the intended physical size.
"""
function plot_closure_comparison(args...; kwargs...)
    _require_glmakie(:plot_closure_comparison)
end

"""
    save_closure_comparison(path, ds, sectors, closures; dpi = 600, kwargs...)

Save a publication PNG, with pixels and PNG resolution metadata matching the
figure's `width_cm` and point-sized text. Default: 17 cm wide at 600 dpi with
10-point text. `kwargs` are forwarded to [`plot_closure_comparison`](@ref).
Returns `path`. Requires `using GLMakie`.
"""
function save_closure_comparison(args...; kwargs...)
    _require_glmakie(:save_closure_comparison)
end
