using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using GenomicBreeding
using GBCore, GBIO, GBModels, GBPlots
using StatsBase, DataFrames
using Dates, ProgressMeter, UnicodePlots
