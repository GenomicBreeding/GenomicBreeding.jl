"""
    writeJLD2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects}; fname::Union{Missing,String} = missing)::String

Save core (`Genomes`, `Phenomes`, and `Trials`) and simulation (`SimulatedEffects`) structs 
as heirarchical data format version 5 (HDF5) - compatible format.
Note that the extension name should be '.jld2'.

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2,p=4); genomes.entries = ["entry_1", "entry_2"]; genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"];

julia> writeJLD2(genomes, fname="test_genomes.jld2")
"test_genomes.jld2"

julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> writeJLD2(phenomes, fname="test_phenomes.jld2")
"test_phenomes.jld2"

julia> trials = Trials(n=1, t=2); trials.entries = ["entry_1"];

julia> writeJLD2(trials, fname="test_trials.jld2")
"test_trials.jld2"

julia> simulated_effects = SimulatedEffects();

julia> writeJLD2(simulated_effects, fname="test_simulated_effects.jld2")
"test_simulated_effects.jld2"
```
"""
function writeJLD2(A::Union{Genomes,Phenomes,Trials,SimulatedEffects}; fname::Union{Missing,String} = missing)::String
    # Check input arguments
    if !checkdims(A)
        throw(DimensionMismatch(string(typeof(A)) * " input is corrupted."))
    end
    if ismissing(fname)
        fname = string("output-", string(typeof(A)), "-", Dates.format(now(), "yyyymmddHHMMSS"), ".jld2")
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if split(basename(fname), ".")[end] != "jld2"
            throw(ArgumentError("The extension name should be `jld2`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Save as a JLD2 binary file
    save(fname, Dict(string(typeof(A)) => A))
    return fname
end


"""
    writedelimited(genomes::Genomes, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Genomes` struct as a string-delimited (default=`\t`) file. 
Each row corresponds to a locus-allele combination.
The first 4 columns correspond to the chromosome, position, all alleles in the locus (delimited by `|`), and the specific allele.
The subsequency columns refer to the samples, pools, entries or genotypes.
Note that the extension name should be '.tsv', '.csv', or '.txt'.

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2,p=4); genomes.entries = ["entry_1", "entry_2"]; genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"];

julia> writedelimited(genomes, fname="test_genomes.tsv")
"test_genomes.tsv"
```
"""
function writedelimited(genomes::Genomes; sep::String = "\t", fname::Union{Missing,String} = missing)::String
    # genomes = Genomes(n=2,p=4); genomes.entries = ["entry_1", "entry_2"]; genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Genomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header line
        header::Array{String,1} = ["chrom", "pos", "all_alleles", "allele"]
        append!(header, genomes.entries)
        header[end] *= "\n"
        write(file, join(header, sep))
        # Rest of the data
        for i = 1:length(genomes.loci_alleles)
            line::Array{String,1} = [string(x) for x in split(genomes.loci_alleles[i], "\t")]
            append!(line, [string(x) for x in genomes.allele_frequencies[:, i]])
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end


"""
    writedelimited(phenomes::Phenomes, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Phenomes` struct as a string-delimited (default=`\t`) file. 
Each row corresponds to a samples, pools, entries or genotypes.
The first 2 columns correspond to the entry and population names.
The subsequency columns refer to the traits containing the phenotype values of each entry.
Note that the extension name should be '.tsv', '.csv', or '.txt'.

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"];

julia> writedelimited(phenomes, fname="test_phenomes.tsv")
"test_phenomes.tsv"
```
"""
function writedelimited(phenomes::Phenomes; sep::String = "\t", fname::Union{Missing,String} = missing)::String
    # phenomes = Phenomes(n=2, t=2); phenomes.entries = ["entry_1", "entry_2"]; phenomes.traits = ["trait_1", "trait_2"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(phenomes)
        throw(DimensionMismatch("Phenomes input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Phenomes-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header line
        header::Array{String,1} = ["entry", "population"]
        append!(header, phenomes.traits)
        header[end] *= "\n"
        write(file, join(header, sep))
        # Rest of the data
        for i in eachindex(phenomes.entries)

            line::Array{String,1} = [phenomes.entries[i], phenomes.populations[i]]
            append!(line, string.(phenomes.phenotypes[i, :]))
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end



"""
    writedelimited(trials::Trials, sep::String = "\t", fname::Union{Missing,String} = missing)::String

Save `Trials` struct as a string-delimited (default=`\t`) file. 
Each row corresponds to a samples, pools, entries or genotypes.
The first 2 columns correspond to the entry and population names.
The subsequency columns refer to the traits containing the phenotype values of each entry.
Note that the extension name should be '.tsv', '.csv', or '.txt'.

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2); trials.years = ["year_1"]; trials.seasons = ["season_1"]; trials.harvests = ["harvest_1"]; trials.sites = ["site_1"]; trials.entries = ["entry_1"]; trials.populations = ["population_1"]; trials.replications = ["replication_1"]; trials.blocks = ["block_1"]; trials.rows = ["row_1"]; trials.cols = ["col_1"]; trials.traits = ["trait_1", "trait_2"];

julia> writedelimited(trials, fname="test_trials.tsv")
"test_trials.tsv"
```
"""
function writedelimited(trials::Trials; sep::String = "\t", fname::Union{Missing,String} = missing)::String
    # trials = Trials(n=1, t=2); trials.years = ["year_1"]; trials.seasons = ["season_1"]; trials.harvests = ["harvest_1"]; trials.sites = ["site_1"]; trials.entries = ["entry_1"]; trials.populations = ["population_1"]; trials.replications = ["replication_1"]; trials.blocks = ["block_1"]; trials.rows = ["row_1"]; trials.cols = ["col_1"]; trials.traits = ["trait_1", "trait_2"]; sep::String = "\t"; fname = missing;
    # Check input arguments
    if !checkdims(trials)
        throw(DimensionMismatch("Trials input is corrupted."))
    end
    if ismissing(fname)
        if sep == "\t"
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".tsv")
        elseif (sep == ",") || (sep == ";")
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".csv")
        else
            fname = string("output-Trials-", Dates.format(now(), "yyyymmddHHMMSS"), ".txt")
        end
    else
        if isfile(fname)
            throw(ErrorException("The file: " * fname * " exists. Please remove or rename the output file."))
        end
        if (split(basename(fname), ".")[end] != "tsv") &&
           (split(basename(fname), ".")[end] != "csv") &&
           (split(basename(fname), ".")[end] != "txt")
            throw(ArgumentError("The extension name should be either `tsv`, `csv` or `txt`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Write into a new text file
    open(fname, "w") do file
        # Header line
        header::Array{String,1} = ["years", "seasons", "harvests", "sites", "entries", "populations", "replications", "blocks", "rows", "cols"]
        append!(header, trials.traits)
        header[end] *= "\n"
        write(file, join(header, sep))
        # Rest of the data
        for i in eachindex(trials.entries)
            line::Array{String,1} = [trials.years[i], trials.seasons[i], trials.harvests[i], trials.sites[i], trials.entries[i], trials.populations[i], trials.replications[i], trials.blocks[i], trials.rows[i], trials.cols[i]]
            append!(line, string.(trials.phenotypes[i, :]))
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    return fname
end