function write(genomes::Genomes, fname::Union{Missing,String} = missing)::String
    # genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false)
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted."))
    end
    if ismissing(fname)
        fname = "output-" * Dates.format(now(), "yyyymmddHHMMSS")  * ".jld2"
    else
        if split(basename(fname), ".")[end] != ".jld2"
            throw(ArgumentError("The extension name should be `.jld2`."))
        end
        if dirname(fname) != ""
            if !isdir(dirname(fname))
                throw(ArgumentError("Directory " *  dirname(fname) * " does not exist."))
            end
        end
    end
    save(fname, Dict("genomes"=>genomes))
    fname
end


function write(genomes::Genomes, sep::String="\t", fname::Union{Missing,String} = missing)::String
    if !checkdims(genomes)
        throw(DimensionMismatch("Genomes input is corrupted."))
    end

    # genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false)
    df::DataFrame = DataFrame(hcat(g.entries, g.populations, g.allele_frequencies), :auto)
    column_names::Array{String,1} = ["entries", "populations"]
    append!(column_names, g.loci_alleles)

end

