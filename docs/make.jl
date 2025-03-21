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
    pages = [
        "Home" => "index.md",
        "Methods Reference" => "references.md",
    ],
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreeding.jl", devbranch = "main")
