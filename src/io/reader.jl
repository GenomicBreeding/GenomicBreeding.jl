function read(type::Type{Genomes}, fname::Union{Missing,String} = missing)::Genomes
    # genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false)
    if !isfile(fname)
        throw(ArgumentError("JLD2 file: " * fname * " does not exist."))
    end
    genomes::Genomes = load(fname)["genomes"]
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes struct from the JLD2 file: " * fname * " is corrupted."))
    end
    return genomes
end

function read(type::Type{Phenomes}, fname::Union{Missing,String} = missing)::Phenomes
    # trials, _effects = simulatetrials(genomes=simulategenomes(verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false)
    # phenomes = Phenomes(n=length(trials.entries), t=length(trials.traits))
    # phenomes.entries = trials.entries
    # phenomes.populations = trials.populations
    # phenomes.traits = trials.traits
    # phenomes.phenotypes = trials.phenotypes
    # checkdims(phenomes)
    if !isfile(fname)
        throw(ArgumentError("JLD2 file: " * fname * " does not exist."))
    end
    genomes::Phenomes = load(fname)["phenomes"]
    if !checkdims(phenomes)
        throw(DimensionMismatch("Phenomes struct from the JLD2 file: " * fname * " is corrupted."))
    end
    return genomes
end
