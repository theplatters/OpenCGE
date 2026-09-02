include("bootstrap.jl")

data = Data("I-O_DE2019_formatiert.csv")
sector = "Vorb.Baustellen-,Bauinstallations-,Ausbauarbeiten"
ci = findfirst(==(sector), data.io.Sektoren)
println("construction index = ", ci)
println("nsectors = ", length(data.consumption_share))
sh = autonomous_shock(data; autonomous_mult = 0.2)
println("autonomous_demand nonzero indices: ", findall(!iszero, sh.autonomous_demand))
println("autonomous_demand[ci] = ", sh.autonomous_demand[ci])
println("consumption_share[ci] = ", data.consumption_share[ci])
println("sum(labor_share) = ", sum(data.labor_share))
