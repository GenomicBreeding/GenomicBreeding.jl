using GenomicBreeding, GBIO, GBModels, GBPlots
using Documenter


DocMeta.setdocmeta!(GenomicBreeding, :DocTestSetup, :(using GenomicBreeding); recursive = true)

makedocs(;
    modules = [GenomicBreeding, GBIO, GBModels, GBPlots],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreeding.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreeding.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Reference" => "references.md",
    ],
    doctest = false,
    checkdocs=:exports,
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreeding.jl", devbranch = "main")
