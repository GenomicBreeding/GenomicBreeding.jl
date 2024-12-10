function read(fname::Union{Missing,String} = missing)::Genomes
    # genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false)
    if !isfile(fname)
        throw(ArgumentError("JLD2 file: " * fname * " does not exist."))
    end
    genomes::Genomes = load(fname)["genomes"]
    if !checkdims(genomes)
        throw(
            DimensionMismatch(
                "Genomes struct from the JLD2 file: " * fname * " is corrupted.",
            ),
        )
    end
    genomes
end
