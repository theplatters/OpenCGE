using CairoMakie
using GraphMakie
using Graphs,SimpleWeightedGraphs
using FileIO
include("../model.jl")

data = CESModel.read_data("I-O_DE2019_formatiert.csv")
img = load(joinpath(pwd(),"data/map_germany.png"))
G = SimpleWeightedDiGraph(data.Ω[1:20,1:20])
edws = [get_weight(G, e.src, e.dst) for e in edges(G)]

f = Figure()
ax = Axis(f[1,1])
hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
image!(ax,rotr90(img))   
graphplot!(ax,G, edge_width= 4 .* log.(1 .+ edws), arrow_size=5 .+ 5 .* edws)
f

