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

include("core/all_structs.jl")
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

export Genomes, Phenomes, Trials, SimulatedEffects, TEBV
export hash, ==, checkdims, dimensions, loci_alleles, loci, plot, slice, filter
export simulategenomes
export simulateeffects, simulategenomiceffects, simulatetrials
export tabularise
export countlevels, @string2formula, trialsmodelsfomulae!, analyse
export writeJLD2, writedelimited, readJLD2, readdelimited


# Precompile
@compile_workload begin
    n = 10
    genomes = simulategenomes(n = n)
    trials, effects = simulatetrials(
        genomes = genomes,
        f_add_dom_epi = [
            0.50 0.25 0.13
            0.90 0.00 0.00
        ],
        n_years = 1,
        n_seasons = 1,
        n_harvests = 1,
        n_sites = 1,
        n_replications = 2,
    )
    phenomes = Phenomes(n = n, t = 2)
    phenomes.entries = trials.entries[1:n]
    phenomes.populations = trials.populations[1:n]
    phenomes.traits = trials.traits
    phenomes.phenotypes = trials.phenotypes[1:n, :]
    phenomes.mask .= true
    tebv = analyse(trials, max_levels = 10)

    readJLD2(Genomes, writeJLD2(genomes))
    readJLD2(Phenomes, writeJLD2(phenomes))
    readJLD2(Trials, writeJLD2(trials))
    readJLD2(SimulatedEffects, writeJLD2(effects[1]))
    readJLD2(TEBV, writeJLD2(tebv))
end

end
