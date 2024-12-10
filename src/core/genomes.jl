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
Genomes([#undef, #undef], [#undef, #undef], [#undef, #undef], Union{Missing, Float64}[missing missing; missing missing], Bool[0 0; 0 0])

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
        new(
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, p),
            fill(missing, n, p),
            fill(false, n, p),
        )
    end
end

"""
    checkdims(genomes::Genomes)::Bool

Check dimension compatibility of the fields of the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2,p=4);

julia> checkdims(genomes)
true

julia> genomes.entries = ["beaking_change"];

julia> checkdims(genomes)
false
```
"""
function checkdims(genomes::Genomes)::Bool
    n, p = size(genomes.allele_frequencies)
    if (n != length(genomes.entries)) ||
       (n != length(genomes.populations)) ||
       (p != length(genomes.loci_alleles)) ||
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
    dimensions(genomes::Genomes)::Tuple{Int64, Int64, Int64}

Count the number of entries, populations, loci, and maximum number of alleles per locus in the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> dimensions(genomes)
(100, 100, 3000, 1000, 4)
```
"""
function dimensions(genomes::Genomes)::Tuple{Int64,Int64,Int64,Int64,Int64}
    n_entries::Int64 = length(genomes.entries)
    n_populations::Int64 = length(genomes.populations)
    n_loci_alleles::Int64 = length(genomes.loci_alleles)
    n_loci::Int64 = 0
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
            n_loci += 1
        else
            if ((chr == locus_ids[1]) && (pos != parse(Int64, locus_ids[2]))) ||
               (chr != locus_ids[1])
                chr = locus_ids[1]
                pos = parse(Int64, locus_ids[2])
                n_alleles = length(split(locus_ids[3], '|'))
                max_n_alleles = max_n_alleles < n_alleles ? n_alleles : max_n_alleles
                n_loci += 1
            end
        end
    end
    (n_entries, n_populations, n_loci_alleles, n_loci, max_n_alleles)
end
