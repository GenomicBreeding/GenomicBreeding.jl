using GenomicBreeding
using Documenter

DocMeta.setdocmeta!(GenomicBreeding, :DocTestSetup, :(using GenomicBreeding); recursive = true)

makedocs(;
    modules = [GenomicBreeding],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreeding.jl",
    format = Documenter.HTML(;
        canonical = "https://jeffersonfparil.github.io/GenomicBreeding.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/jeffersonfparil/GenomicBreeding.jl", devbranch = "main")
