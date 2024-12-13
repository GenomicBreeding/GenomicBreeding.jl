module GenomicBreeding

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
# Plots.backend(:plotly)
Plots.backend(:gr)

include("core/genomes.jl")
include("core/phenomes.jl")
include("core/trials.jl")
include("simulation/simulate_effects.jl")
include("simulation/simulate_genomes.jl")
include("simulation/simulate_trials.jl")
include("simulation/simulate_mating.jl")
include("precompile.jl")

export Genomes, Phenomes, Trials, SimulatedEffects
export checkdims, dimensions, loci_alleles, loci, plot, slice
export simulategenomes
export simulateeffects, simulategenomiceffects, simulatetrials
export tabularise

end
