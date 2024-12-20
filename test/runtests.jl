using GenomicBreeding, Test, Documenter

Documenter.doctest(GenomicBreeding)

# TODO: Add more unit tests in addition to the doctests
@testset "GenomicBreeding.jl" begin
    genomes::Genomes = simulategenomes(; verbose = false)
    @test isa(genomes, Genomes)
end
