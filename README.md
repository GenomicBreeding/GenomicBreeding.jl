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