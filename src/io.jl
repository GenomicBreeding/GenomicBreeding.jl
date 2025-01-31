"""
    mutable struct GBInput
        fname_geno::String
        fname_pheno::String
        bulk_cv::Bool
        populations::Union{Nothing, String, Vector{String}}
        traits::Union{Nothing, String, Vector{String}}
        models::Any
        n_folds::Int64
        n_replications::Int64
        maf::Float64
        mtv::Float64
        verbose::Bool
    end

Input struct

- `fname_geno`: genotype file (see file format guide: **TODO:** {URL HERE})
- `fname_pheno`: phenotype file (see file format guide: **TODO:** {URL HERE})
- `bulk_cv`: perform cross-validation across all populations, i.e. disregards population grouping (Default = false)
- `populations`: include only these populations (Default = nothing which means include all populations)
- `traits`: include only these traits (Default = nothing which means include all traits)
- `models`: include these model functions (Default = [ridge, bayesa]; see models list: **TODO:** {URL HERE})
- `n_folds`: number of k partitions for k-fold cross-validation (Default = 5)
- `n_replications`: number of replications for repeated k-fold cross-validation (Default = 5)
- `maf`: minimum allele frequency (Default = 0.05)
- `mtv`: minimum trait variance (Default = 1e-7)
- `verbose`: show messages (Default = true)
"""
mutable struct GBInput
    fname_geno::String
    fname_pheno::String
    bulk_cv::Bool
    populations::Union{Nothing, String, Vector{String}}
    traits::Union{Nothing, String, Vector{String}}
    models::Any
    n_folds::Int64
    n_replications::Int64
    maf::Float64
    mtv::Float64
    verbose::Bool
    function GBInput(;fname_geno::String,
        fname_pheno::String,
        bulk_cv::Bool = false,
        populations::Union{Nothing, String, Vector{String}} = nothing,
        traits::Union{Nothing, String, Vector{String}} = nothing,
        models::Any = [ridge, bayesa],
        n_folds::Int64 = 5,
        n_replications::Int64 = 5,
        maf::Float64 = 0.05,
        mtv::Float64 = 1e-7,
        verbose::Bool = true,
    )
        new(
            fname_geno,
            fname_pheno,
            bulk_cv,
            populations,
            traits,
            models,
            n_folds,
            n_replications,
            maf,
            mtv,
            verbose,
        )
    end
end

"""
    load(input::GBInput)::Tuple{Genomes, Phenomes}

Load, merge and filter genotype and phenotype data

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding)
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = writedelimited(genomes, fname="test-geno.tsv");

julia> fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv");

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], verbose=false);

julia> genomes, phenomes = load(input);

julia> length(unique(genomes.populations)) == length(unique(phenomes.populations)) == 2
true

julia> length(phenomes.traits) == 1
true
"""
function load(input::GBInput)::Tuple{Genomes, Phenomes}
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = writedelimited(genomes, fname="test-geno.tsv")
    # # fname_geno = writejld2(genomes, fname="test-geno.jld2")
    # # fname_geno = writevcf(genomes, fname="test-geno.vcf") ### Will have unknown population groupings
    # fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv")
    # # fname_pheno = writejld2(phenomes, fname="test-pheno.jld2")
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # Parse input
    fname_geno = input.fname_geno
    fname_pheno = input.fname_pheno
    bulk_cv = input.bulk_cv
    populations = input.populations
    traits = input.traits
    # models = input.models
    n_folds = input.n_folds
    # n_replications = input.n_replications
    maf = input.maf
    mtv = input.mtv
    verbose = input.verbose
    # Load genomes and phenomes
    genomes = try
        readdelimited(Genomes, fname=fname_geno)
    catch
        try
            readjld2(Genomes, fname=fname_geno)
        catch
            try
                readvcf(fname=fname_geno)
                @warn "The vcf file: `" * fname_geno * "` lacks population grouping information. This will likely result in conflicts in merging with th phenotype data."
            catch
                throw(ArgumentError("Unrecognised genotype file format: `" * fname_geno * "`.\n" * 
                "Please refer to the file format guide: **TODO:** {URL HERE}"))
            end
        end
    end
    phenomes = try
        readdelimited(Phenomes, fname=fname_pheno)
    catch
        try
            readjld2(Phenomes, fname=fname_pheno)
        catch
            throw(ArgumentError("Unrecognised phenotype file format: `" * fname_pheno * "`.\n" * 
            "Please refer to the file format guide: **TODO:** {URL HERE}"))
        end
    end
    # Merge the genomes and phenomes
    genomes, phenomes = merge(genomes, phenomes)
    # Slice by population and traits
    populations = if isnothing(populations)
        unique(sort(genomes.populations))
    else
        populations
    end
    traits = if isnothing(traits)
        sort(phenomes.traits)
    else
        traits
    end
    # Check population names
    all_populations = unique(sort(genomes.populations))
    unrecognised_populations = []
    for pop in populations
        if !(pop ∈ all_populations)
            push!(unrecognised_populations, pop)
        end
    end
    if length(unrecognised_populations) > 0
        throw(ArgumentError("Unrecognised population/s:\n\t‣ " * join(unrecognised_populations, "\n\t‣ ")))
    end
    # Check trait names
    unrecognised_traits = []
    for trait in traits
        if !(trait ∈ phenomes.traits)
            push!(unrecognised_traits, trait)
        end
    end
    if length(unrecognised_traits) > 0
        throw(ArgumentError("Unrecognised trait/s:\n\t‣ " * join(unrecognised_traits, "\n\t‣ ")))
    end
    # Check if we have enough entries per fold
    if !bulk_cv
        for population in populations
            n = sum(genomes.populations .== population)
            if (n/n_folds) < 10
                throw(ArgumentError("The numer of entries in population: `" * population * "` results in less than 10 entries per " *
                string(n_folds) *"-fold cross-validation. " * 
                "Please consider reducing the number of folds or setting `cv_bulk = true`."))
            end
        end
    end
    idx_traits = findall([sum(traits .== trait) > 0 for trait in phenomes.traits])
    idx_entries = findall([sum(populations .== pop) > 0 for pop in genomes.populations])
    n = length(idx_entries)
    if (n / n_folds) < 10
        throw(ArgumentError("The number of entries left after filtering by populations results in less than 10 entries per " *
        string(n_folds) *"-fold cross-validation. " * 
        "Please consider reducing the number of folds or including more populations."))
    end
    # Slice phenomes to retain only the requested population/s and/or trait/s
    if length(populations) < length(all_populations)
        phenomes = slice(phenomes, idx_entries=idx_entries)
    end
    if length(traits) < length(phenomes.traits)
        phenomes = slice(phenomes, idx_traits=idx_traits)
    end
    # Filter phenomes by mtv, i.e. retain traits with variance
    fixed_traits = []
    for trait in phenomes.traits
        # trait = phenomes.traits[1]
        if var(phenomes.phenotypes[:, phenomes.traits .== trait]) < mtv
            push!(fixed_traits, trait)
        end
    end
    if length(fixed_traits) > 0
        @warn "Excluding the following fixed trait/s:\n\t‣ " * join(fixed_traits, "\n\t‣ ")
        idx_traits = findall([sum(fixed_traits .== trait) == 0 for trait in phenomes.traits])
        if length(idx_traits) == 0
            throw(ArgumentError("All requested traits are fixed, i.e.\n\t‣ " * join(phenomes.traits, "\n\t‣ ") * "\n. On the following requested populations:\n\t‣ " * join(populations, "\n\t‣ ")))
        end
        phenomes = slice(phenomes, idx_traits=idx_traits)
    end
    # Slice the genomes to retain only the requested populations
    if length(populations) < length(all_populations)
        genomes = slice(genomes, idx_entries=idx_entries)
    end
    # Filter genomes by maf (minimum allele frequency)
    genomes = filter(genomes, maf)
    # Show dimensions of the input genomes and phenomes after merging and filterings
    if verbose
        println("Genomes: ")
        @show dimensions(genomes)
        println("Phenomes: ")
        @show dimensions(phenomes)
    end
    # Output
    if genomes.entries != phenomes.entries
        throw(ErrorException("Error loading the genomes and phenomes data from:\n\t‣ Genotype file: " * fname_geno * "\n\t‣ Phenotype file: " * fname_pheno))
    end
    (genomes, phenomes)
end

