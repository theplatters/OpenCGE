using JET

@report_opt data = Data("I-O_DE2019_formatiert.csv")
@report_opt Model(data, shocks, ces_options)


@report_opt solve(Model(data, shocks, ces_options))