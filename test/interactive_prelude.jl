using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using GenomicBreeding
using StatsBase, DataFrames
using Dates, ProgressMeter, UnicodePlots
