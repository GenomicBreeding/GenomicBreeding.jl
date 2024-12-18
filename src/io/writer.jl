function writeJLD2(
    A::Union{Genomes,Phenomes,Trials,SimulatedEffects};
    fname::Union{Missing,String} = missing,
)::String
    # A = Genomes()
    # A = Phenomes()
    # A = Trials()
    # A = SimulatedEffects()
    # Check input arguments
    if !checkdims(A)
        throw(DimensionMismatch(string(typeof(A)) * " input is corrupted."))
    end
    if ismissing(fname)
        fname = string(
            "output-",
            string(typeof(A)),
            "-",
            Dates.format(now(), "yyyymmddHHMMSS"),
            ".jld2",
        )
    else
        if isfile(fname)
            throw(
                ErrorException(
                    "The file: " *
                    fname *
                    "exists. Please remove or rename the output file.",
                ),
            )
        end
        if split(basename(fname), ".")[end] != ".jld2"
            throw(ArgumentError("The extension name should be `.jld2`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " * dirname(fname) * " does not exist."))
            end
        end
    end
    # Save as a JLD2 binary file
    save(fname, Dict(string(typeof(A)) => A))
    fname
end


function writedelimited(
    genomes::Genomes,
    sep::String = "\t",
    fname::Union{Missing,String} = missing,
)::String
    # genomes = simulategenomes(n=100, l=1_000); sep::String = "\t"; fname::Union{Missing,String} = missing;
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
            throw(
                ErrorException(
                    "The file: " *
                    fname *
                    "exists. Please remove or rename the output file.",
                ),
            )
        end
        if (split(basename(fname), ".")[end] != ".tsv") ||
           (split(basename(fname), ".")[end] != ".csv") ||
           (split(basename(fname), ".")[end] != ".txt")
            throw(ArgumentError("The extension name should be `.tsv`, `.csv` or `.txt`."))
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
            line::Array{String,1} =
                [string(x) for x in split(genomes.loci_alleles[i], "\t")]
            append!(line, [string(x) for x in genomes.allele_frequencies[:, i]])
            line[end] *= "\n"
            write(file, join(line, sep))
        end
    end
    fname
end
