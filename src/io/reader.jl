using DataFrames, CSV

ret = open("data.dat", "r") do datafile
    [parse.(Float64, split(line)) for line in eachline(datafile)]
end
