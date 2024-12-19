module GenomicBreeding

using PrecompileTools: @compile_workload

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

include("core/genomes.jl")
include("core/phenomes.jl")
include("core/trials.jl")

include("simulation/simulate_effects.jl")
include("simulation/simulate_genomes.jl")
include("simulation/simulate_trials.jl")
include("simulation/simulate_mating.jl")

include("models/tebv.jl")
include("models/ols.jl")

include("io/reader.jl")
include("io/writer.jl")

# include("precompile.jl")

export Genomes, Phenomes, Trials, SimulatedEffects, TEBV
export hash, ==, checkdims, dimensions, loci_alleles, loci, plot, slice
export simulategenomes
export simulateeffects, simulategenomiceffects, simulatetrials
export tabularise
export countlevels, @string2formula, trialsmodelsfomulae!, analyse
export writeJLD2, writedelimited


# Precompile
@compile_workload begin
    g = Genomes()
    p = Phenomes()
    t = Trials()
    e = SimulatedEffects()
    checkdims(g)
    checkdims(p)
    checkdims(t)
    checkdims(e)

    genomes = simulategenomes(; n = 2, l = 2, n_chroms = 1, verbose = false)
    dimensions(genomes)
    simulateeffects()
    simulategenomiceffects(; genomes = genomes)
    trials, vector_of_effects = simulatetrials(genomes = genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false)
    tebv = analyse(trials, max_levels=10, verbose=false);

end

end
