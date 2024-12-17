# GenomicBreeding

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jeffersonfparil.github.io/GenomicBreeding.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeffersonfparil.github.io/GenomicBreeding.jl/dev/)
[![Build Status](https://github.com/jeffersonfparil/GenomicBreeding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jeffersonfparil/GenomicBreeding.jl/actions)
[![Coverage](https://codecov.io/gh/jeffersonfparil/GenomicBreeding.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jeffersonfparil/GenomicBreeding.jl)

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
    user="jeffersonfparil",
    authors=["jeffersonparil@gmail.com"],
    dir="./SimQuantGen.jl",
    julia=v"1.11",
    plugins=[
        License(; name="GPL-3.0+", path=nothing, destination="LICENSE.md"),
        CompatHelper(),
        Codecov(),    
        GitHubActions(;
        osx=false,
        windows=false,
        ),
        Documenter{GitHubActions}(),
        Git(;
            ignore=["*.jl.*.cov",
                    "*.jl.cov",
                    "*.jl.mem",
                    "*.code-workspace",
                    ".DS_Store",
                    "docs/build/",
                    "Manifest.toml",
                    "tmp/",
                    "*.svg"],
            ssh=true
        ),
    ],
)
t("SimQuantGen.jl")
```