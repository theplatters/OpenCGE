using MAT

function loadInData()
    file = matopen("data/simulationData.mat")

    data = read(file, "data")

    close(file)

    return data
end

function loadStfp()
    file = matopen("data/stfp.mat")

    stfp = read(file,"stfp")

    close(file)

    return stfp
end