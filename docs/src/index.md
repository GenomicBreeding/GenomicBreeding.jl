# GenomicBreeding.jl

Documentation for the `GenomicBreeding` module.

```@contents
```

## Overview

The `GenomicBreeding` module provides a comprehensive suite of tools for genomic prediction, genome-wide association studies (GWAS), and data handling in genomic breeding. It integrates functionalities from `GBCore`, `GBIO`, `GBModels`, and `GBPlots` to offer efficient and scalable solutions for genetic data analysis and visualisation.

## Installation

We designed [GenomicBreeding.jl](https://github.com/GenomicBreeding/GenomicBreeding.jl) to work on an HPC running Linux (the various components, i.e. [GBCore.jl](https://github.com/GenomicBreeding/GBCore.jl), [GBIO.jl](https://github.com/GenomicBreeding/GBIO.jl), [GBModels.jl](https://github.com/GenomicBreeding/GBModels.jl), and [GBPlots.jl](https://github.com/GenomicBreeding/GBPlots.jl) work on a single Linux PC too).

Currently, we require that you install [Julia](https://julialang.org/) on your home directory in your HPC cluster via:

```shell
curl -fsSL https://install.julialang.org | sh
type -a julia
```

Currently, [GBModels.jl](https://github.com/GenomicBreeding/GBModels.jl) is dependent on [R](https://www.r-project.org/) and the package [BGLR](https://github.com/gdlc/BGLR-R) for Bayes A, Bayes B and Bayes C models. Because of this we require that [R](https://www.r-project.org/) and [BGLR](https://github.com/gdlc/BGLR-R) be installed. To help with this, you may install all the requirements via [Conda](https://www.anaconda.com/docs/getting-started/miniconda/main) using the environment file: [`GenomicBreeding_conda.yml`](https://github.com/GenomicBreeding/GenomicBreeding.jl/blob/main/GenomicBreeding_conda.yml). We aim to have a pure Julia implementation of Bayesian models using [Turing.jl](https://turinglang.org/) in the near future (we just need to speed-up the models a bit).


Install the [GenomicBreeding.jl](https://github.com/GenomicBreeding/GenomicBreeding.jl) library in Julia:

```julia
using Pkg
Pkg.add("https://github.com/GenomicBreeding/GenomicBreeding.jl")
```

Feel free to install the [GenomicBreeding.jl components](https://github.com/GenomicBreeding) as well as various other useful libraries:

```julia
using Pkg
GB_components = [
    "https://github.com/GenomicBreeding/GBCore.jl",
    "https://github.com/GenomicBreeding/GBIO.jl",
    "https://github.com/GenomicBreeding/GBModels.jl",
    "https://github.com/GenomicBreeding/GBPlots.jl",
]
for P in GB_components
    Pkg.add(url=P)
end
Pkg.add(["StatsBase", "MixedModels", "MultivariateStats", "UnicodePlots", "ColorSchemes", "CairoMakie"])
```

## Quickstart

- [Genomic prediction](#genomic-prediction)
- [Genome-wide association study](#genome-wide-association-study)
- [Genotype data filtering and imputation](#genotype-data-filtering-and-imputation)
- [Phenotype data analyses](#phenotype-data-analyses)
- [Cluster analyses](#cluster-analyses)
- [Plotting](#plotting)
- [Mating simulations](#mating-simulations)

### Genomic prediction

#### Example 1: using simulated data

Here's a simple example using simulated data to perform replicated k-fold cross validation:

```julia
# It is always a good idea to keep all your packages updated
using Pkg; Pkg.update()
# Load GenomicBreeding
using GenomicBreeding
# Load plotting and GP model functions
import GenomicBreeding: plot, lasso, bayesa
# Simulate genotype and phenotype data
genomes = simulategenomes(n=300, l=1_000, verbose=true)
genomes.populations[1:100] .= "pop1"; genomes.populations[101:200] .= "pop2"; genomes.populations[201:300] .= "pop3" # simulate multiple populations
trials, _ = simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=true);
phenomes = extractphenomes(trials)
fname_geno = writedelimited(genomes, fname="test-geno.tsv")
fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv")
# Input struct documentation
@doc GBInput
# Setup the input struct
input = GBInput(
    fname_geno=fname_geno, 
    fname_pheno=fname_pheno,
    anaysis=cv, # set analysis to use the `cv` function for replicated k-fold cross-validation
    models = [lasso, bayesa],
    n_folds=2, 
    n_replications=2, 
    SLURM_job_name="testGB",
    SLURM_account_name="dbiopast1",
    SLURM_cpus_per_task=2, 
    SLURM_mem_G=5, 
    SLURM_time_limit_dd_hhmmss="00-00:15:00",
    SLURM_max_array_jobs_running=10,
    verbose=true
)
# Preliminary look at the genotype and phenotype data
plot(input=input, format="png", plot_size=(700, 500))
# Documentation of the main user interface function (take note of the currently available analyses)
@doc submitslurmarrayjobs
# Submit the Slurm array jobs
# Note that here you will be prompted to enter YES to proceed.
outdir = submitslurmarrayjobs(input)
# Monitor the Slurm jobs
run(`sh -c 'squeue -u $USER'`)
run(`sh -c 'ls -lhtr slurm-*_*.out'`)
run(`sh -c 'cat slurm-*_*.out'`)
run(`sh -c 'tail slurm-*_*.out'`)
run(`sh -c 'grep -i "err" slurm-*_*.out'`)
run(`sh -c 'grep -i "err" slurm-*_*.out | cut -d: -f1 | sort | uniq'`)
readdir(outdir)
# Once the array jobs have finishes or at least a couple of jobs have finished, run below.
# Rerun as you often as wish to update the plots.
# You may exit Julia and just run the plotting function below after correctly defining input::GBInput above.
plot(input=input, format="png", plot_size=(700, 500), skip_genomes=true, skip_phenomes=true, overwrite=true)
```

#### Example 2: using empirical data

```julia
# TODO
# with GP per se on unphenotyped set
```

### Genome-wide association study

```julia
# TODO
```

### Genotype data filtering and imputation

```julia
# TODO
```

### Phenotype data analyses

```julia
# TODO
```

### Cluster analyses

```julia
# TODO
```

### Plotting

```julia
# TODO
```

### Mating simulations

```julia
# TODO
```

## License

The `GenomicBreeding` module is licensed under the GPLv3 License. See the [LICENSE.md](https://github.com/GenomicBreeding/GenomicBreeding.jl/blob/main/LICENSE.md) file for more details.

## [Index](@id main-index)

```@index
Pages = ["references.md"]
```
