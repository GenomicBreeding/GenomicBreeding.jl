# Manual

```@contents
Pages = ["manual.md"]
Depth = 3
```

## Main

- [`GBInput`](https://genomicbreeding.github.io/GenomicBreeding.jl/dev/references/#GenomicBreeding.GBInput): Struct for input data.
- `checkinputs`: Check input data.
- `loadgenomesphenomes`: Load genomic and phenomic data.
- `loadcvs`: Load cross-validation data.
- `loadfits`: Load model fits.
- `prepareinputs`: Prepare input data.
- `prepareoutprefixandoutdir`: Prepare output prefix and directory.
- `cv`: Perform cross-validation.
- `fit`: Fit models.
- `predict`: Make predictions.
- `gwas`: Perform genome-wide association studies.
- `submitslurmarrayjobs`: Submit SLURM array jobs.
- `plot`: General plotting function.

## Core Data Structures

- `AbstractGB`: Abstract type for genomic breeding data.
- `Genomes`: Struct for genomic data.
- `Phenomes`: Struct for phenomic data.
- `Trials`: Struct for trials data.
- `SimulatedEffects`: Struct for simulated genetic effects.
- `TEBV`: Struct for trial-estimated breeding values.
- `Fit`: Struct for genotype-to-phenotype models.
- `CV`: Struct for cross-validation of genotype-to-phenotype models.

## General Functions

- `clone`: Clone a struct.
- `hash`: Compute the hash of a struct.
- `==`: Check equality of structs.
- `checkdims`: Check dimensions of genomic data.
- `dimensions`: Get dimensions of genomic data.
- `loci_alleles`: Get loci alleles.
- `loci`: Get loci.
- `plot`: Plot genomic data.
- `slice`: Slice genomic data.
- `filter`: Filter genomic data.
- `tabularise`: Convert genomic data to a tabular format.
- `summarise`: Summarise genomic data.

## Simulation Functions

- `simulategenomes`: Simulate genomic data.
- `simulateeffects`: Simulate genetic effects.
- `simulategenomiceffects`: Simulate genomic effects.
- `simulatetrials`: Simulate trials data.

## Phenotype Analysis Functions

- `countlevels`: Count levels in phenomic data.
- `@string2formula`: Convert string to formula.
- `trialsmodelsfomulae!`: Generate trial models formulae.
- `analyse`: Analyse phenomic data.
- `extractphenomes`: Extract phenomic data.
- `@stringevaluation`: Evaluate string expressions.
- `addcompositetrait`: Add composite trait to phenomic data.

## Input/Output Functions (GBIO)

- `levenshteindistance`: Compute Levenshtein distance for fuzzy matching.
- `isfuzzymatch`: Check for fuzzy matches.
- `readjld2`: Read JLD2 files.
- `readdelimited`: Read delimited text files.
- `readvcf`: Read VCF files.
- `writejld2`: Write JLD2 files.
- `writedelimited`: Write delimited text files.
- `writevcf`: Write VCF files.

## Genomic Prediction Functions (GBModels)

- `metrics`: Compute metrics for model evaluation.
- `extractxyetc`: Extract features and targets for modeling.
- `predict`: Make predictions using trained models.
- `grmsimple`: Compute genomic relationship matrix (simple).
- `grmploidyaware`: Compute genomic relationship matrix (ploidy-aware).
- `gwasprep`: Prepare data for GWAS.
- `gwasols`: Perform GWAS using ordinary least squares.
- `gwaslmm`: Perform GWAS using linear mixed models.
- `loglikreml`: Compute log-likelihood for REML.
- `gwasreml`: Perform GWAS using REML.
- `square`, `invoneplus`, `log10epsdivlog10eps`, `mult`, `addnorm`, `raise`: Utility functions for transformations.
- `transform1`, `transform2`, `epistasisfeatures`, `@string2operations`, `reconstitutefeatures`: Functions for feature engineering.
- `bglr`, `bayesian`: Bayesian genomic prediction models.
- `turing_bayesG`, `turing_bayesGs`, `turing_bayesGπ`, `turing_bayesGπs`, `turing_bayesL`, `turing_bayesLs`, `turing_bayesLπ`, `turing_bayesLπs`, `turing_bayesT`, `turing_bayesTπ`, `turing_bayesG_logit`: Turing Bayesian models.
- `ols`, `ridge`, `lasso`, `bayesa`, `bayesb`, `bayesc`: Linear and Bayesian regression models.
- `validate`: Validate models.
- `cvmultithread!`, `cvbulk`: Cross-validation functions.
- `cvperpopulation`, `cvpairwisepopulation`, `cvleaveonepopulationout`: Cross-validation strategies.

## Plotting Functions (GBPlots)

- `PlotsGB`: General plotting functions.
- `DistributionPlots`: Distribution plots.
- `ViolinPlots`: Violin plots.
- `CorHeatPlots`: Correlation heatmaps.
- `TreePlots`: Tree plots.
- `BarPlots`: Bar plots.
- `BoxPlots`: Box plots.
- `PCBiPlots`: Principal component biplots.
- `checkdims`: Check dimensions for plotting.
- `labeltofname`: Convert labels to filenames.
- `saveplots`: Save plots to files.
