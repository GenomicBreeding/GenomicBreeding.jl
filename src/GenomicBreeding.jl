module GenomicBreeding

using GBCore
export AbstractGB, Genomes, Phenomes, Trials, SimulatedEffects, TEBV, Fit, CV
export clone, hash, ==
export checkdims, dimensions, loci_alleles, loci, plot, slice, filter, tabularise, summarise
export simulategenomes, simulateeffects, simulategenomiceffects, simulatetrials
export countlevels, @string2formula, trialsmodelsfomulae!, analyse, extractphenomes
export @stringevaluation, addcompositetrait

using GBIO
export levenshteindistance, isfuzzymatch
export readjld2, readdelimited, readvcf
export writejld2, writedelimited, writevcf

using GBModels
export metrics
export extractxyetc, predict
export grmsimple, grmploidyaware
export gwasprep, gwasols, gwaslmm, loglikreml, gwasreml
export square, invoneplus, log10epsdivlog10eps, mult, addnorm, raise
export transform1, transform2, epistasisfeatures, @string2operations, reconstitutefeatures
export bglr, bayesian
export turing_bayesG, turing_bayesGs, turing_bayesGπ, turing_bayesGπs
export turing_bayesL, turing_bayesLs, turing_bayesLπ, turing_bayesLπs
export turing_bayesT, turing_bayesTπ
export turing_bayesG_logit
export ols, ridge, lasso, bayesa, bayesb, bayesc
export validate, cvmultithread!, cvbulk
export cvperpopulation, cvpairwisepopulation, cvleaveonepopulationout

using GBPlots
export PlotsGB, DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots
export checkdims, labeltofname, saveplots

using StatsBase, DataFrames, Dates
using ProgressMeter, UnicodePlots


include("io.jl")
include("genomic_prediction.jl")
include("gwas.jl")
include("slurm.jl")
include("plot.jl")

export GBInput, hash, ==, clone, checkinputs, loadgenomesphenomes, loadcvs, loadfits, prepareinputs, prepareoutprefixandoutdir
export cv, fit, predict, gwas
export submitslurmarrayjobs
export plot



end
