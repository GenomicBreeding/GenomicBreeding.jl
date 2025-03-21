# GenomicBreeding

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions)

## Installation

We designed [GenomicBreeding.jl](https://github.com/GenomicBreeding/GenomicBreeding.jl) to work on an HPC running Linux (the various components, i.e. [GBCore.jl](https://github.com/GenomicBreeding/GBCore.jl), [GBIO.jl](https://github.com/GenomicBreeding/GBIO.jl), [GBModels.jl](https://github.com/GenomicBreeding/GBModels.jl), and [GBPlots.jl](https://github.com/GenomicBreeding/GBPlots.jl) work on a single Linux PC too).

Currently, we require that you install Julia on your home directory in you HPC cluster via:

```shell
curl -fsSL https://install.julialang.org | sh
type -a julia
```

Currently, [GBModels.jl](https://github.com/GenomicBreeding/GBModels.jl) is dependent on [R](https://www.r-project.org/) and the package [BGLR](https://github.com/gdlc/BGLR-R) for Bayes A, Bayes B and Bayes C models. Because of this we require that [R](https://www.r-project.org/) and [BGLR](https://github.com/gdlc/BGLR-R) be installed. To help with this, you may install all the requirements via [Conda](https://www.anaconda.com/docs/getting-started/miniconda/main) using the environment file: [`GenomicBreeding_conda.yml`](GenomicBreeding_conda.yml). We aim to have a pure Julia implementation of Bayesian models using [Turing.jl](https://turinglang.org/) in the near future (we just need to speed-up the models a bit).


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

### 1. Example 1: simulated data

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

### 2. Example 2: test data

```julia

```

## File formats

### Variant call format 

See [VCFv4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and [VCFv4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf) for details in the format specifications.

Note that given that we work with autopolyploids, and pools (e.g. half-sib families) more often, the VCF parser prioritise the `AF` (allele frequency) field, followed by the `AD` (allele depth) field, and finally the `GT` (genotype) field. If your VCF file needs depth-based filtering, you may opt for having the `AD` field instead of the `AF` field, as well as include the `DP` (depth) field. A way to generate a VCF file with `AD` and `DP` fields is via: 

```shell
time \
bcftools mpileup \
    -BI \
    -a AD,DP \
    -d 4000000 \
    -f ${REF} \
    -T ${LOCI} \
    -b ${BAMS_LIST} \
    -Ov > ${PREFIX}-MPILEUP.vcf
```

where the following flags refer to:

- `-B`: disable probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.
- `-I`: do not perform INDEL calling
- `-a AD,DP`: comma-separated list of FORMAT and INFO tags to output, i.e. Total allelic depth (Number=R,Type=Integer) (AD) and Number of high-quality bases (Number=1,Type=Integer) (DP)
- `-Ou`: output uncompressed BCF, because we're piping into bcftools call
- `-d 4000000`: at each position, read maximally 4 million reads per input file
- `-f ...`: reference genome in fasta format
- `-T ...`: regions or loci specified in a tab-delimited file. The columns of the tab-delimited file can contain either positions (two-column format: CHROM, POS) or intervals (three-column format: CHROM, BEG, END), but not both. Positions are 1-based and inclusive.
- `-b ...`: list of input alignment files, one file per line

### Allele frequency table

This is a simple human-readable string-delimited text file (tab-delimited by default). The smallest minimal valid allele frequency table is as follows:

| chrom | pos | all_alleles | allele | entry_1 |
| :---- | :-- | :---------- | :----- | :------ |
| chr_1 | 123 | A|T         | A      | 0.5     |

Each column must be sorted exactly as this: starting with "**chrom**" for chromosome or scaffold name, "**pos**" for the position in bases, "**all_alleles**" for a string of reference and alternative alleles separated by pipe/s (`|`), "**allele**" for exactly one of the alleles in the previous column, and finally subsequency column names refer to the names of the entries. Values under the 5th and subsequent columns are assumed to be allele frequencies, i.e. ranges from `0.0` to `1.0`. Missing values may be encoded as any of the following:

- ""
- "missing"
- "NA"
- "na"
- "N/A"
- "n/a"

An additional header may be included, for example:

| chrom | pos | all_alleles | allele | entry_1 | entry_2 | entry_3 |
| chrom | pos | all_alleles | allele | pop_ABC | pop_DEF | pop_XYZ |
| :---- | :-- | :---------- | :----- | :------ | :------ | :------ |
| chr_1 | 123 | A\|T        | A      | 0.51    | 0.25    | 0.25    |

where the first header is the same as the detailed above, while the second header replaces the entry names with the corresponding population or group the entries belong to.

### JLD2

This is a compressed binary format containing Julia structs (e.g. genomes, and phenomes struct). It is a subset of the scientific data format [HDF5](https://www.hdfgroup.org/solutions/hdf5/).

### Trials table

Similar to the allele frequency table format, this is a simple human-readable string-delimited text file (tab-delimited by default). The smallest minimal valid trials table is as follows:

| years | seasons           | harvests          | sites     | entries | populations | replications | blocks | rows | cols | trait_1 |
| :---- | :---------------- | :---------------- | :-------- | :------ | :---------- | :----------- | :----- | :--- | :--- | :------ |
| 2023  | 2023_Early_Spring | FLIGHT-2023-09-05 | LOC1_TRT1 | entry_1 | pop_1       | rep_1        | row1   | row1 | col1 | 3.1416  |

The first 10 columns must be named (or as close as possible) as the following in any order: 

- "**years**"
- "**seasons**"
- "**harvests**"
- "**sites**"
- "**entries**"
- "**populations**"
- "**replications**"
- "**blocks**"
- "**rows**"
- "**cols**"

The subsequent columns (column 11 and so on) refer to the name of the traits. Values under the 11th and subsequent columns are assumed to be numeric. Missing values may be encoded as any of the following:

- ""
- "missing"
- "NA"
- "na"
- "N/A"
- "n/a"


### Phenomes table

Again, similar to the allele frequency table format, this is a simple human-readable string-delimited text file (tab-delimited by default). The smallest minimal valid phenomes table is as follows:

| entries | populations | trait_1 |
| :------ | :---------- | :------ |
| entry_1 | pop_XYZ     | 0.5772  |

Each column must be sorted exactly as this: starting with "**entries**" for the names of the entries, genotypes, families or pools, "**populations**" for the corresponding population names or group names. Subsequent column names (column 3 and so on) refer to the trait names. Values under the 3rd and subsequent columns are assumed to be numeric. Missing values may be encoded as any of the following:

- ""
- "missing"
- "NA"
- "na"
- "N/A"
- "n/a"

## Dev stuff:

<details>
<summary>Details</summary>

### REPL prelude

```shell
julia --threads 8,1 --load test/interactive_prelude.jl
```

### Format and test

```shell
time julia test/cli_tester.jl
```

### Initialise a new package

```julia
using PkgTemplates
t = Template(;
    user="GenomicBreeding",
    authors=["jeffersonparil@gmail.com"],
    dir="./",
    julia=v"1.11",
    plugins=[
        License(; name="GPL-3.0+", path=nothing, destination="LICENSE.md"),
        CompatHelper(),
        GitHubActions(;
        osx=false,
        windows=false,
        ),
        Documenter{GitHubActions}(),
        Git(;
            ignore=[
                "*.code-workspace",
                "*.jl.*.cov",
                "*.jl.cov",
                "*.jl.mem",
                ".DS_Store",
                "/docs/Manifest.toml",
                "/docs/build/",
                "Manifest.toml",
                "docs/build/",
                "tmp/",
                "*.svg",
                "*.jld2",
                "*.tsv",
                "*.csv",
                "*.txt"
            ],
            ssh=true
        ),
    ],
)
t("GBPlots.jl")
```

Install slurm:

```shell
sudo apt install slurmd slurmctld -y
sudo chmod 777 /etc/slurm
sudo cat << EOF > /etc/slurm/slurm.conf
# slurm.conf file generated by configurator.html.
# Put this file on all nodes of your cluster.
# See the slurm.conf man page for more information.
#
ClusterName=localcluster
SlurmctldHost=localhost
MpiDefault=none
ProctrackType=proctrack/linuxproc
ReturnToService=2
SlurmctldPidFile=/var/run/slurmctld.pid
SlurmctldPort=6817
SlurmdPidFile=/var/run/slurmd.pid
SlurmdPort=6818
SlurmdSpoolDir=/var/lib/slurm/slurmd
SlurmUser=slurm
StateSaveLocation=/var/lib/slurm/slurmctld
SwitchType=switch/none
TaskPlugin=task/none
#
# TIMERS
InactiveLimit=0
KillWait=30
MinJobAge=300
SlurmctldTimeout=120
SlurmdTimeout=300
Waittime=0
# SCHEDULING
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core
#
#AccountingStoragePort=
AccountingStorageType=accounting_storage/none
JobCompType=jobcomp/none
JobAcctGatherFrequency=30
JobAcctGatherType=jobacct_gather/none
SlurmctldDebug=info
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdDebug=info
SlurmdLogFile=/var/log/slurm/slurmd.log
#
# COMPUTE NODES
NodeName=localhost CPUs=1 RealMemory=500 State=UNKNOWN
PartitionName=LocalQ Nodes=ALL Default=YES MaxTime=INFINITE State=UP
EOF
sudo chmod 755 /etc/slurm/
sudo systemctl start slurmctld
sudo systemctl start slurmd
sudo scontrol update nodename=localhost state=idle
sinfo
sudo cat /var/log/slurm/slurmd.log
sudo cat /var/log/slurm/slurmctld.log
```

Install Lmod:

```shell
sudo apt install lua5.4 liblua5.4-dev lmod -y
sudo apt install tcl-dev -y
wget https://sourceforge.net/projects/lmod/files/Lmod-8.7.tar.bz2
tar xfvj Lmod-8.7.tar.bz2
rm Lmod-8.7.tar.bz2
cd Lmod-8.7/
./configure --prefix=$HOME --with-fastTCLInterp=no
sudo make install
echo 'export PATH=$HOME/lmod/8.7/libexec:$PATH' >> ~/.bashrc
echo 'source $HOME/lmod/8.7/init/bash' >> ~/.bashrc
echo 'export LMOD_CMD=$HOME/lmod/8.7/libexec/lmod' >> ~/.bashrc
echo 'export MODULEPATH="/etc/lmod/modules/"' >> ~/.bashrc
```

Sample module file (`/etc/lmod/modules/R.lua`):

```shell
sudo chmod -R 777 /etc/lmod/modules/
sudo cat << EOF > /etc/lmod/modules/R.lua
help([[
...
]])
whatis("Version: 4.1.2")
whatis("R statistical computing environment")
prepend_path("LD_LIBRARY_PATH","/usr/local/lib/R/site-library/")
prepend_path("LIBRARY_PATH","\$HOME/R/x86_64-pc-linux-gnu-library/4.3")
prepend_path("PATH","/usr/bin")
EOF
sudo chmod -R 755 /etc/lmod/modules/
```

Test

```shell
module avail R
module add R
```
</details>
