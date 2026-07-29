# ── Data loader (extracted from test_standalone.jl) ──

function read_io_table(filename)
    raw = read(filename, String)
    raw = replace(raw, "\r\n" => "\n")
    lines = split(raw, "\n")
    header_line = lines[1] * "\n" * lines[2]
    header = split(header_line, ';')
    header[1] = replace(header[1], "\"" => "")
    n_cols = length(header)
    data_lines = filter(!isempty, lines[3:end])
    sectors = String[]
    data_matrix = zeros(length(data_lines), n_cols - 1)
    for (i, row) in enumerate(data_lines)
        parts = split(row, ';')
        push!(sectors, strip(parts[1]))
        for j in 2:min(n_cols, length(parts))
            val_str = replace(parts[j], "," => ".")
            val_str = strip(val_str)
            if val_str in ["-", "x", ""]
                data_matrix[i, j-1] = 0.0
            else
                data_matrix[i, j-1] = parse(Float64, val_str)
            end
        end
    end
    return sectors, header[2:end], data_matrix
end

function load_model_data()
    sectors, colnames, io_matrix = read_io_table("data/I-O_DE2019_formatiert.csv")
    colnames = strip.(colnames)
    sectors = strip.(sectors)

    N = 71
    bws_row = findfirst(==("Bruttowertschöpfung"), sectors)
    prodwert_row = findfirst(==("Produktionswert"), sectors)
    arbeit_row = findfirst(==("Arbeitnehmerentgelt im Inland"), sectors)
    gesamt_verwendung_col = findfirst(==("Gesamte Verwendung von Gütern"), colnames)
    konsum_col = findfirst(==("Konsumausgaben der privaten Haushalte im Inland"), colnames)
    exporte_col = findfirst(==("Exporte"), colnames)

    Ω = io_matrix[1:N, 1:N]
    Ω = Ω ./ sum(Ω, dims=2)
    Ω[isnan.(Ω)] .= 0.0

    grossy = io_matrix[1:N, gesamt_verwendung_col]
    value_added = io_matrix[bws_row, 1:N]
    labor_comp = io_matrix[arbeit_row, 1:N]

    factor_share = value_added ./ grossy
    factor_share = clamp.(factor_share, 0.01, 0.99)

    consumption = zeros(N)
    for col in konsum_col:exporte_col
        consumption .+= io_matrix[1:N, col]
    end
    consumption_share_gross_output = consumption ./ grossy

    consumption_share_raw = (I - diagm(1 .- factor_share) * Ω)' * grossy
    consumption_share_raw[consumption_share_raw .< 0] .= 0.0
    consumption_share = consumption_share_raw ./ sum(consumption_share_raw)

    λ = inv(I - diagm(1 .- factor_share) * Ω)' * consumption_share
    labor_share = λ .* factor_share

    # Standard shock
    target_sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"
    sector_idx = findfirst(==(target_sector), sectors[1:N])
    demand_shock = ones(N)
    supply_shock = ones(N)
    demand_shock[sector_idx] = 1.8097957577943152

    return (
        Ω=Ω, consumption_share=consumption_share, factor_share=factor_share,
        supply_shock=supply_shock, demand_shock=demand_shock,
        labor_share=labor_share, λ=λ, N=N, grossy=grossy,
        sectors=sectors, konsum_col=konsum_col, exporte_col=exporte_col,
        value_added=value_added, labor_comp=labor_comp,
    )
end