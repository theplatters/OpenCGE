function load_impulses(filename)
  filedir = joinpath(pwd(), "data/", filename)
  impulses = CSV.read(filedir, DataFrames.DataFrame, delim=",", decimal='.', missingstring=["-", "x"]) #read in from csv
  select!(impulses, Not(names(impulses)[1])) # remove the first column
  impulses ./ 1_000_000
end
