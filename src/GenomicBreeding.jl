module GenomicBreeding

using Random
using StatsBase
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
# Plots.backend(:plotly)
Plots.backend(:gr)

include("core/genomes.jl")
include("core/phenomes.jl")
include("core/trials.jl")
include("simulation/simulate_effects.jl")
include("simulation/simulate_genomes.jl")
include("simulation/simulate_trials.jl")
include("precompile.jl")

export Genomes, Phenomes, Trials, SimulatedEffects
export checkdims, dimensions, loci_alleles, loci, slice
export simulategenomes
export simulateeffects, simulategenomiceffects, simulatetrials
export tabularise, plot

end
