# GenomicBreeding

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions)

## Installation

We designed [GenomicBreeding.jl](https://github.com/GenomicBreeding/GenomicBreeding.jl) to work on an HPC running Linux (the various components, i.e. [GBCore.jl](https://github.com/GenomicBreeding/GBCore.jl), [GBIO.jl](https://github.com/GenomicBreeding/GBIO.jl), [GBModels.jl](https://github.com/GenomicBreeding/GBModels.jl), and [GBPlots.jl](https://github.com/GenomicBreeding/GBPlots.jl) work on a single Linux PC too).

Currently, we require that you install Julia on your home directory in you HPC cluster via:

```shell
curl -fsSL https://install.julialang.org | sh
type -a julia
```

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

Here's a simple example using simulated data:

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
    models = [lasso, bayesa],
    n_folds=2, 
    n_replications=2, 
    SLURM_job_name="testGB",
    SLURM_account_name="dbiopast1",
    SLURM_cpus_per_task=2, 
    SLURM_mem_G=5, 
    SLURM_time_limit_dd_hhmmss="00-00:15:00",
    SLURM_max_array_jobs_running=10,
    SLURM_module_load_R_version_name="R",
    verbose=true
)
# Preliminary look at the genotype and phenotype data
plot(input=input, format="png", plot_size=(700, 500))
# Perform replicated k-fold cross-validation
outdir = submitslurmarrayjobs(input=input, analysis=assess)
# Monitor the Slurm jobs
run(`sh -c 'squeue -u $USER'`)
run(`sh -c 'ls -lhtr slurm-*_*.out'`)
run(`sh -c 'cat slurm-*_*.out'`)
run(`sh -c 'tail slurm-*_*.out'`)
run(`sh -c 'grep -i "err" slurm-*_*.out'`)
run(`sh -c 'grep -i "err" slurm-*_*.out | cut -d: -f1 | sort | uniq'`)
readdir(outdir)
# Once the array jobs have finishes or at least a couple of jobs have finished, run below and rerun as you wish to update the plots:
plot(input=input, format="png", plot_size=(700, 500), skip_genomes=true, skip_phenomes=true, overwrite=true)
```

### 2. Example 2: test data


## Dev stuff:

<details>
<summary>Details</summary>

### REPL prelude

```shell
julia --threads 8,1 --load test/prelude.jl
```

### Format and test

```shell
time julia test/cli_tester.jl
```

### Docstring conventions

- Structs and main functions with title description, etc including Examples with doctests
- Methods, i.e. functions with the same names but different input types follow the usual Julia docstring pattern, i.e. the function signature, then some description, then details including parameter description, and finally examples with doctests

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