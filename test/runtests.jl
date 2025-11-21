using GenomicBreeding
using Test
using Documenter, DocumenterTools

try
    Documenter.doctest(GenomicBreeding)
catch
    DocumenterTools.generate()
    Pkg.activate("docs/")
    Documenter.doctest(GenomicBreeding)
end

@testset "GenomicBreeding.jl" begin
    @test 1 == 1
end
