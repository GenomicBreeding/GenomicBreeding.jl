# GenomicBreeding

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GenomicBreeding.github.io/GenomicBreeding.jl/dev/)
[![Build Status](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GenomicBreeding/GenomicBreeding.jl/actions)

## Dev stuff:

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

```

Install Lmod:

```shell
wget https://sourceforge.net/projects/lmod/files/lua-5.1.4.9.tar.bz2
tar xf lua-5.1.4.9.tar.bz2
cd lua-5.1.4.9
./configure --prefix=/opt/apps/lua/5.1.4.9
make 
sudo make install
cd /opt/apps/lua
sudo ln -s 5.1.4.9 lua
sudo ln -s /opt/apps/lua/lua/bin/lua /usr/local/bin
rm lua-5.1.4.9.tar.bz2
sudo apt install lua5.4 liblua5.4-dev lmod -y
sudo apt -y install tcl-dev

wget https://sourceforge.net/projects/lmod/files/Lmod-8.7.tar.bz2
tar xfvj Lmod-8.7.tar.bz2
rm Lmod-8.7.tar.bz2
cd Lmod-8.7/
./configure --prefix=$HOME
sudo make install
export PATH=$HOME/lmod/8.7/libexec:$PATH
source $HOME/lmod/8.7/init/bash
export LMOD_CMD=$HOME/lmod/8.7/libexec/lmod

export MODULEPATH="/etc/lmod/modules/"
```

Sample module file (`/etc/lmod/modules/applications/R.lua`):

```lua
help([[
For detailed instructions, go to:
   https://...

]])
whatis("Version: 4.3")
whatis("R statistical computing environment")
setenv("R_HOME", "/usr/bin/R")
setenv("R_LIBS_USER", "$HOME/R/x86_64-pc-linux-gnu-library/4.3")
prepend_path("PATH", "/usr/bin")

alias("R", "/usr/bin/R")
```

Test

```shell
module avail
module add applications/R
```
