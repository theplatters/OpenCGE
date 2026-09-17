
# Plotting entry points.
#
# The plotting implementation lives in the GLMakie package extension
# (`ext/BeyondHultenGLMakieExt.jl`), which Julia loads automatically when
# `GLMakie` and `BeyondHulten` are both loaded in the same session. The
# stubs below keep the exported plotting names defined for headless use and
# raise an instructive error telling the user how to enable plotting.

function _require_glmakie(fname::Symbol)
	error("`BeyondHulten.$(fname)` requires the GLMakie extension, which is not loaded. " *
		"Run `using GLMakie` together with `using BeyondHulten` " *
		"(install it first with `import Pkg; Pkg.add(\"GLMakie\")` if necessary).")
end

"""
	axis_change_in_level!(fig, data, impulses; options)

Left panel of `panel`: bar chart of demand shocks with error bars.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function axis_change_in_level!(args...; kwargs...)
	_require_glmakie(:axis_change_in_level!)
end

"""
	axis_change_in_price!(fig, data, impulse; options)

Right panel of `panel`: price/quantity scatter.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function axis_change_in_price!(args...; kwargs...)
	_require_glmakie(:axis_change_in_price!)
end

"""
	get_color(data, shocks)

Highlight vector for the shocked sectors in `axis_change_in_price!`.
Plot helper implemented by the GLMakie extension; requires `using GLMakie`.
"""
function get_color(args...; kwargs...)
	_require_glmakie(:get_color)
end

"""
	panel(data, impulses; options, name = "panel")

Two-panel impulse-response figure (bar chart + price/quantity scatter).
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function panel(args...; kwargs...)
	_require_glmakie(:panel)
end

"""
	diff_lambda(data, impulses; options, name = "diff_lambda_imp")

Stacked CGE vs. Leontief sectoral decomposition figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function diff_lambda(args...; kwargs...)
	_require_glmakie(:diff_lambda)
end

"""
	effect_of_different_elasticities(shocks, data, gdp_effect_simple; labor_slack_function, name)

Real-GDP elasticity-gradient figure for several elasticity levels.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function effect_of_different_elasticities(args...; kwargs...)
	_require_glmakie(:effect_of_different_elasticities)
end

"""
	comparison_between_labor_slacks(data, shocks, gdp_effect_simple, title)

GDP figure comparing labour-slack specifications.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function comparison_between_labor_slacks(args...; kwargs...)
	_require_glmakie(:comparison_between_labor_slacks)
end

"""
	labor_slack_gradient(data, impulse)

Real-GDP figure along the labour-slack interpolation.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function labor_slack_gradient(args...; kwargs...)
	_require_glmakie(:labor_slack_gradient)
end

"""
	plot_real_gdp_gradient(results; title, cd, leontief, initial, ylims)

Four-panel real-GDP elasticity-gradient figure.
Implemented by the GLMakie extension; requires `using GLMakie`.
"""
function plot_real_gdp_gradient(args...; kwargs...)
	_require_glmakie(:plot_real_gdp_gradient)
end
