using Pkg
Pkg.activate(".")
using GenomicBreeding

using Random
using StatsBase
using StatsModels
using MixedModels
using Distributions
using LinearAlgebra
using ProgressMeter
using DataFrames
using FileIO
using JLD2
using CSV
using Dates
using UnicodePlots
using Plots
using StatsPlots
using Metida
# Plots.backend(:plotly)
Plots.backend(:gr)