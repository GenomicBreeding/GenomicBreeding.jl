using GenomicBreeding, Test, Documenter

Documenter.doctest(GenomicBreeding)

@testset "GenomicBreeding.jl" begin
    genomes::Genomes = simulategenomes(verbose = false)
    @test isa(genomes, Genomes)
end
