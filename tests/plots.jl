using Plots
using BeyondHulten
using CSV, DataFrames
using Plotly

plotly()

f(x, y) = min(x/0.55,y/0.45)

x = range(1, 2.5, length=100)
y = range(1, 2.5, length=50)
z = @. f(x', y)
a = Plots.contour(x, y, z, title ="Leontief")

e = 0.2

f(x,y) = (0.55 * x^((e-1)/e) + 0.45 * y^((e-1)/e))^(e/(e-1))
x = range(0, 2.5, length=100)
y = range(0, 2.5, length=50)
z = @. f(x', y)
b = contour(x, y, z, title = "CES ϵ = 0.2")

e = 0.7

f(x,y) = (0.55 * x^((e-1)/e) + 0.45 * y^((e-1)/e))^(e/(e-1))
x = range(0, 2.5, length=100)
y = range(0, 2.5, length=50)
z = @. f(x', y)
c = contour(x, y, z, title = "CES ϵ = 0.7")

f(x,y) = (x^0.55 * y^0.45)
x = range(0, 2.5, length=100)
y = range(0, 2.5, length=50)
z = @. f(x', y)
d = contour(x, y, z, title = "Cobb-Douglas")

plot(a,b,c,d, size=(800,600))
Plots.savefig("plots/production_functions.png")

GDP_nominal = CSV.read("data/demand_shock_nominal.csv",DataFrame)
GDP_no_realloc_nominal = CSV.read("data/demand_shock_nominal_no_realloc.csv",DataFrame)
GDP = CSV.read("data/demand_shock.csv",DataFrame)
GDP_no_realloc = CSV.read("data/demand_shock_no_realloc.csv",DataFrame)

min_shock = 1.0
max_shock = 1.8
shock_count = 50
lay = @layout [a b; c d]

p1 = Plots.scatter(title ="GDP")
p2 = Plots.scatter(title="GDP no realloc")
p3 = Plots.scatter(title="GDP nominal")
p4 = Plots.scatter(title = "GDP no realloc nominal")
for sector in names(select(GDP,Not(:DemandShockAmount)))
    lab = replace(sector, "Dienstleistungen" => "DL")
    Plots.scatter!(p1, LinRange(min_shock, max_shock, shock_count), GDP[!,sector], label=lab, legend=false)
    Plots.scatter!(p2,LinRange(min_shock, max_shock, shock_count), GDP_no_realloc[!,sector], label=lab, legend=false)
    Plots.scatter!(p3, LinRange(min_shock, max_shock, shock_count), GDP_nominal[!,sector], label=lab, legend=false)
    Plots.scatter!(p4,LinRange(min_shock, max_shock, shock_count), GDP_no_realloc_nominal[!,sector], label=lab, legend=false)
end

p1

Plots.plot(p1, p2, p3, p4, layout=lay, legend = false, size=(1980, 1080))

Plots.savefig("plots/shocks.html")


data = Data("I-O_DE2019_formatiert.csv")

data.io 