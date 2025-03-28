module GenomicBreeding

using GenomicBreedingCore
export AbstractGB, Genomes, Phenomes, Trials, SimulatedEffects, TEBV, Fit, CV, GRM
export clone, hash, ==
export checkdims, dimensions, loci_alleles, loci, distances, plot, tabularise, summarise
export slice, sparsities, filter, filterbysparsity, filterbymaf, filterbypca, filterbysnplist
export simulatechromstruct,
    simulateposandalleles, simulatepopgroups, simulateldblocks, simulateperpopμΣ, simulateallelefreqs!
export simulategenomes, simulateeffects, simulategenomiceffects, simulatetrials
export histallelefreqs, simulatemating
export countlevels, @string2formula, trialsmodelsfomulae!, analyse, extractphenomes
export @stringevaluation, addcompositetrait
export maskmissing!, divideintomockscaffolds, estimateld, estimatedistances, knni, knnioptim, impute
export inflatediagonals!, grmsimple, grmploidyaware

using GenomicBreedingIO
export levenshteindistance, isfuzzymatch
export readjld2, writejld2
export readdelimited, writedelimited
export vcfcountloci,
    vcfchunkify,
    vcfextractentriesandformats,
    vcfextractinfo,
    vcfinstantiateoutput,
    vcfparsecoordinates,
    vcfextractallelefreqs!
export readvcf, writevcf


using GenomicBreedingModels
export metrics
export extractxyetc, predict
export gwasprep, gwasols, gwaslmm, loglikreml, gwasreml
export square, invoneplus, log10epsdivlog10eps, mult, addnorm, raise
export transform1, transform2, epistasisfeatures, @string2operations, reconstitutefeatures
export bglr, bayesian
export turing_bayesG, turing_bayesGs, turing_bayesGπ, turing_bayesGπs
export turing_bayesL, turing_bayesLs, turing_bayesLπ, turing_bayesLπs
export turing_bayesT, turing_bayesTπ
export turing_bayesG_logit
export mlp
export ols, ridge, lasso, bayesa, bayesb, bayesc
export validate, cvmultithread!, cvbulk
export cvperpopulation, cvpairwisepopulation, cvleaveonepopulationout

using GenomicBreedingPlots
export PlotsGB, DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots, BarPlots, BoxPlots, PCBiPlots
export checkdims, labeltofname, saveplots
export plot

using StatsBase, DataFrames, Dates
using ProgressMeter, UnicodePlots

include("io.jl")
include("genomic_prediction.jl")
include("gwas.jl")
include("slurm.jl")
include("plot.jl")

export GBInput, clone, hash, ==
export checkinputs, loadgenomesphenomes, loadcvs, loadfits, prepareinputs, prepareoutprefixandoutdir
export cv, fit, predict, gwas
export submitslurmarrayjobs
export plot

end
