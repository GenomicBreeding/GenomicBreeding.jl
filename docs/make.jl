using Pkg
Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
Pkg.add(url = "https://github.com/GenomicBreeding/GBIO.jl")
Pkg.add(url = "https://github.com/GenomicBreeding/GBModels.jl")
Pkg.add(url = "https://github.com/GenomicBreeding/GBPlots.jl")
Pkg.develop(url = "https://github.com/GenomicBreeding/GenomicBreeding.jl")
Pkg.add("Documenter")
using GenomicBreeding
using Documenter


DocMeta.setdocmeta!(GenomicBreeding, :DocTestSetup, :(using GenomicBreeding); recursive = true)

makedocs(;
    modules = [GenomicBreeding],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreeding.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreeding.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreeding.jl", devbranch = "main")
