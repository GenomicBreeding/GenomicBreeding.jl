"""
# Genomes struct 

Containes unique entries and loci_alleles where allele frequencies can have missing values

## Constructor

```julia
Genomes(; n::Int64 = 1, p::Int64 = 2)
```

## Fields
- `entries`: names of the `n` entries or samples
- `populations`: name/s of the population/s each entries or samples belong to
- `loci_alleles`: names of the `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles) including the chromsome or scaffold name, position, all alleles, and current allele separated by tabs
- `allele_frequencies`: `n x p` matrix of allele frequencies between 0 and 1 which can have missing values
- `mask`: `n x p` matrix of boolean mask for selective analyses and slicing

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2, p=2)
Genomes(["", ""], ["", ""], ["", ""], Union{Missing, Float64}[missing missing; missing missing], Bool[0 0; 0 0])

julia> fieldnames(Genomes)
(:entries, :populations, :loci_alleles, :allele_frequencies, :mask)

julia> genomes.entries = ["entry_1", "entry_2"];

julia> genomes.populations = ["pop_1", "pop_1"];

julia> genomes.loci_alleles = ["chr1_12345_A|T_A", "chr2_678910_C|D_D"];

julia> genomes.allele_frequencies = [0.5 0.25; 0.9 missing];

julia> genomes.mask = [true true; true false];

julia> genomes
Genomes(["entry_1", "entry_2"], ["pop_1", "pop_1"], ["chr1_12345_A|T_A", "chr2_678910_C|D_D"], Union{Missing, Float64}[0.5 0.25; 0.9 missing], Bool[1 1; 1 0])
```
"""
mutable struct Genomes
    entries::Array{String,1}
    populations::Array{String,1}
    loci_alleles::Array{String,1}
    allele_frequencies::Array{Union{Float64,Missing},2}
    mask::Array{Bool,2}
    function Genomes(; n::Int64 = 1, p::Int64 = 2)
        new(fill("", n), fill("", n), fill("", p), fill(missing, n, p), fill(false, n, p))
    end
end

"""
    checkdims(genomes::Genomes)::Bool

Check dimension compatibility of the fields of the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2,p=4);

julia> checkdims(genomes)
false

julia> genomes.entries = ["entry_1", "entry_2"];

julia> genomes.loci_alleles = ["locus_1", "locus_2", "locus_3", "locus_4"];

julia> checkdims(genomes)
true
```
"""
function checkdims(genomes::Genomes)::Bool
    n, p = size(genomes.allele_frequencies)
    if (n != length(genomes.entries)) ||
       (n != length(unique(genomes.entries))) ||
       (n != length(genomes.populations)) ||
       (p != length(genomes.loci_alleles)) ||
       (p != length(unique(genomes.loci_alleles))) ||
       ((n, p) != size(genomes.mask))
        return false
    end
    if !isa(genomes.entries, Array{String,1}) ||
       !isa(genomes.populations, Array{String,1}) ||
       !isa(genomes.loci_alleles, Array{String,1}) ||
       !isa(genomes.allele_frequencies, Array{Union{Float64,Missing},2}) ||
       !isa(genomes.mask, Array{Bool,2})
        return false
    end
    return true
end

"""
    dimensions(genomes::Genomes)::Dict{String, Int64}

Count the number of entries, populations, loci-alleles combination, loci, and maximum number of alleles per locus in the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> dimensions(genomes)
Dict{String, Int64} with 6 entries:
  "n_entries"      => 100
  "n_chr"          => 7
  "n_loci"         => 1000
  "n_loci_alleles" => 3000
  "n_populations"  => 1
  "max_n_alleles"  => 4
```
"""
function dimensions(genomes::Genomes)::Dict{String,Int64}
    n_entries::Int64 = length(unique(genomes.entries))
    n_populations::Int64 = length(unique(genomes.populations))
    n_loci_alleles::Int64 = length(genomes.loci_alleles)
    n_loci::Int64 = 0
    n_chr::Int64 = 0
    max_n_alleles::Int64 = 0
    chr::String = ""
    pos::Int64 = 0
    for locus in genomes.loci_alleles
        # locus = genomes.loci_alleles[1]
        locus_ids::Array{String,1} = split(locus, '\t')
        if n_loci == 0
            chr = locus_ids[1]
            pos = parse(Int64, locus_ids[2])
            max_n_alleles = length(split(locus_ids[3], '|'))
            n_chr += 1
            n_loci += 1
        else
            if ((chr == locus_ids[1]) && (pos != parse(Int64, locus_ids[2]))) ||
               (chr != locus_ids[1])
                if chr != locus_ids[1]
                    n_chr += 1
                end
                chr = locus_ids[1]
                pos = parse(Int64, locus_ids[2])
                n_alleles = length(split(locus_ids[3], '|'))
                max_n_alleles = max_n_alleles < n_alleles ? n_alleles : max_n_alleles
                n_loci += 1

            end
        end
    end
    Dict(
        "n_entries" => n_entries,
        "n_populations" => n_populations,
        "n_loci_alleles" => n_loci_alleles,
        "n_chr" => n_chr,
        "n_loci" => n_loci,
        "max_n_alleles" => max_n_alleles,
    )
end

"""
loci_alleles(genomes::Genomes)::Tuple{Array{String,1},Array{Int64,1},Array{String,1},}

Extract chromosomes, positions, and alleles across loci-allele combinations

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> chromsomes, positions, alleles = loci_alleles(genomes);

julia> length(chromsomes), length(positions), length(alleles)
(3000, 3000, 3000)
```
"""
function loci_alleles(
    genomes::Genomes,
)::Tuple{Array{String,1},Array{Int64,1},Array{String,1}}
    chromosomes::Array{String,1} = []
    positions::Array{Int64,1} = []
    alleles::Array{String,1} = []
    for locus in genomes.loci_alleles
        # locus = genomes.loci_alleles[1]
        locus_ids::Array{String,1} = split(locus, '\t')
        push!(chromosomes, locus_ids[1])
        push!(positions, parse(Int64, locus_ids[2]))
        push!(alleles, locus_ids[4])
    end
    (chromosomes, positions, alleles)
end


"""
loci(genomes::Genomes)::Tuple{Array{String,1},Array{Int64,1}}

Extract chromosome names, positions, start and end indexes of each locus across loci

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> chromsomes, positions, loci_ini_idx, loci_fin_idx = loci(genomes);

julia> length(chromsomes), length(positions), length(loci_ini_idx), length(loci_fin_idx)
(1000, 1000, 1000, 1000)
```
"""
function loci(
    genomes::Genomes,
)::Tuple{Array{String,1},Array{Int64,1},Array{Int64,1},Array{Int64,1}}
    chromosomes::Array{String,1} = []
    positions::Array{Int64,1} = []
    loci_ini_idx::Array{Int64,1} = []
    loci_fin_idx::Array{Int64,1} = []
    idx::Int64 = 0
    for locus in genomes.loci_alleles
        # locus = genomes.loci_alleles[1]
        idx += 1
        locus_ids::Array{String,1} = split(locus, '\t')
        if length(chromosomes) == 0
            push!(chromosomes, locus_ids[1])
            push!(positions, parse(Int64, locus_ids[2]))
            push!(loci_ini_idx, idx)
        else
            if (
                (chromosomes[end] == locus_ids[1]) &&
                (positions[end] != parse(Int64, locus_ids[2]))
            ) || (chromosomes[end] != locus_ids[1])
                push!(chromosomes, locus_ids[1])
                push!(positions, parse(Int64, locus_ids[2]))
                push!(loci_ini_idx, idx)
                push!(loci_fin_idx, idx - 1)
            end
        end
    end
    if loci_fin_idx[end] < length(genomes.loci_alleles)
        push!(loci_fin_idx, length(genomes.loci_alleles))
    end
    (chromosomes, positions, loci_ini_idx, loci_fin_idx)
end


"""
    slice(genomes::Genomes)::Tuple{Int64, Int64, Int64}

Count the number of entries, populations, loci, and maximum number of alleles per locus in the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> sliced_genomes = slice(genomes, idx_entries=collect(1:10); idx_loci_alleles=collect(1:300));

julia> dimensions(sliced_genomes)
Dict{String, Int64} with 6 entries:
  "n_entries"      => 10
  "n_chr"          => 1
  "n_loci"         => 100
  "n_loci_alleles" => 300
  "n_populations"  => 1
  "max_n_alleles"  => 4
```
"""
function slice(
    genomes::Genomes;
    idx_entries::Array{Int64,1},
    idx_loci_alleles::Array{Int64,1},
)::Genomes
    # genomes::Genomes = simulategenomes(); idx_entries::Array{Int64,1}=sample(1:100, 10); idx_loci_alleles::Array{Int64,1}=sample(1:10_000, 1000);
    genomes_dims::Dict{String,Int64} = dimensions(genomes)
    n_entries::Int64 = genomes_dims["n_entries"]
    n_loci_alleles::Int64 = genomes_dims["n_loci_alleles"]
    if (minimum(idx_entries) < 1) || (maximum(idx_entries) > n_entries)
        throw(ArgumentError("We accept `idx_entries` from 1 to `n_entries` of `genomes`."))
    end
    if (minimum(idx_loci_alleles) < 1) || (maximum(idx_loci_alleles) > n_loci_alleles)
        throw(
            ArgumentError(
                "We accept `idx_loci_alleles` from 1 to `n_loci_alleles` of `genomes`.",
            ),
        )
    end
    sort!(idx_entries)
    sort!(idx_loci_alleles)
    unique!(idx_entries)
    unique!(idx_loci_alleles)
    n, p = length(idx_entries), length(idx_loci_alleles)
    sliced_genomes::Genomes = Genomes(n = n, p = p)
    for (i1, i2) in enumerate(idx_entries)
        sliced_genomes.entries[i1] = genomes.entries[i2]
        sliced_genomes.populations[i1] = genomes.populations[i2]
        for (j1, j2) in enumerate(idx_loci_alleles)
            if i1 == 1
                sliced_genomes.loci_alleles[j1] = genomes.loci_alleles[j2]
            end
            sliced_genomes.allele_frequencies[i1, j1] = genomes.allele_frequencies[i2, j2]
            sliced_genomes.mask[i1, j1] = genomes.mask[i2, j2]
        end
    end
    ### Check dimensions
    if !checkdims(sliced_genomes)
        throw(DimensionMismatch("Error slicing the genome."))
    end
    # Output
    return (sliced_genomes)
end
