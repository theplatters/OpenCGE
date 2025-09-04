#=============================================================================
Imports
===============================================================================#
using BeyondHulten
using GLMakie
using StatsBase
using ThreadsX
using ProgressMeter
using Statistics

const cd_elasticities = CESElasticities(0.99, 0.99, 0.99)
const cd_options = CES(cd_elasticities, model -> model.data.labor_share)
const cd_options_ls = CES(cd_elasticities, BeyondHulten.full_labor_slack)
const ces_elasticities = CESElasticities(0.001, 0.5, 0.9)
const ces_options = CES(ces_elasticities, model -> model.data.labor_share, false)
const ces_options_ls = CES(ces_elasticities, BeyondHulten.full_labor_slack, false)
const ces_options_ls_alt = CES(ces_elasticities, BeyondHulten.full_labor_slack_alt, false)
const leontief = Leontief()
const cd_options_ls_empirical = CES(cd_elasticities, BeyondHulten.empirical_labor_slack)
const ces_options_ls_empirical = CES(ces_elasticities, BeyondHulten.empirical_labor_slack)

function (@main)(args)
	data = Data("I-O_DE2019_formatiert.csv")
	impulses = load_impulses("impulses.csv")

	shocks = impulse_shock(data, impulses)

	gdp_effect_simple = 1 + sum(shocks.demand_shock_raw) ./ sum(data.io[findfirst(==("Bruttowertschöpfung"), data.io.Sektoren), 2:72])

	panel(data, impulses, options = ces_options)
	panel(data, impulses, options = ces_options_ls, name = "panel_ls")
	panel(data, impulses, options = ces_options_ls_alt, name = "panel_ls_alt")
	panel(data, impulses, options = ces_options_ls_empirical, name = "panel_ls_empirical")

	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = (model -> model.data.labor_share), name = "impulse")
	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.full_labor_slack, name = "impulse_ls")
	effect_of_different_elasticities(shocks, data, gdp_effect_simple, labor_slack_function = BeyondHulten.empirical_labor_slack, name = "impulse_ls_empirical")

	comparison_between_labor_slacks(data, shocks, gdp_effect_simple, "")

	diff_lambda(data, impulses, options = ces_options)
	diff_lambda(data, impulses, options = ces_options_ls, name = "diff_lambda_imp_ls")
	diff_lambda(data, impulses, options = ces_options_ls_alt, name = "diff_lambda_imp_ls_alt")
	diff_lambda(data, impulses, options = ces_options_ls_empirical, name = "diff_lambda_imp_ls_empirical")
	labor_slack_gradient(data, impulses)
	open("data/sector_names.txt", "w") do file
		for name in names(impulses)[2:72]
			println(file, name)
		end
	end
end

