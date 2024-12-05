"""
# Genomes struct 

Containes unique entries and loci_alleles where allele frequencies can have missing values

## Constructor

```julia
Genomes(; n::Int = 1, p::Int = 2)
```

## Fields
- `entries`: names of the `n` entries or samples
- `loci_alleles`: names of the `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles) including the chromsome or scaffold name, position, all alleles, and current allele separated by tabs
- `allele_frequencies`: `n x p` matrix of allele frequencies between 0 and 1 which can have missing values
- `mask`: `n x p` matrix of boolean mask for selective analyses and slicing

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2, p=2)
Genomes([#undef, #undef], [#undef, #undef], Union{Missing, Real}[#undef #undef; #undef #undef], Bool[0 0; 0 0])

julia> fieldnames(Genomes)
(:entries, :loci_alleles, :allele_frequencies, :mask)

julia> genomes.entries = ["entry_1", "entry_2"];

julia> genomes.loci_alleles = ["chr1_12345_A|T_A", "chr2_678910_C|D_D"];

julia> genomes.allele_frequencies = [0.5 0.25; 0.9 missing];

julia> genomes.mask = [true true; true false];

julia> genomes
Genomes(["entry_1", "entry_2"], ["chr1_12345_A|T_A", "chr2_678910_C|D_D"], Union{Missing, Real}[0.5 0.25; 0.9 missing], Bool[1 1; 1 0])
```
"""
mutable struct Genomes
    entries::Array{String,1}
    loci_alleles::Array{String,1}
    allele_frequencies::Array{Union{Real,Missing},2}
    mask::Array{Bool,2}
    function Genomes(; n::Int = 1, p::Int = 2)
        new(
            Array{String,1}(undef, n),
            Array{String,1}(undef, p),
            Array{Real,2}(undef, n, p),
            fill(false, n, p),
        )
    end
end

"""
    checkdims(g::Genomes)::Bool

Check dimension compatibility of the fields of the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> g = Genomes(n=2,p=4);

julia> checkdims(g)
true

julia> g.entries = ["beaking_change"];

julia> checkdims(g)
false
```
"""
function checkdims(g::Genomes)::Bool
    n, p = size(g.allele_frequencies)
    if (n != length(g.entries)) || (p != length(g.loci_alleles)) || ((n, p) != size(g.mask))
        return false
    end
    if !isa(g.entries, Array{String,1}) ||
       !isa(g.loci_alleles, Array{String,1}) ||
       !isa(g.allele_frequencies, Array{Union{Real,Missing},2}) ||
       !isa(g.mask, Array{Bool,2})
        return false
    end
    return true
end

"""
    dimensions(g::Genomes)::Tuple{Int, Int, Int}

Count the number of entries, loci, and maximum number of alleles per locus in the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> g = simulategenomes(n=100, l=1_000, n_alleles=4);

julia> dimensions(g)
(100, 3000, 1000, 4)
```
"""
function dimensions(g::Genomes)::Tuple{Int,Int,Int,Int}
    n_entries::Int = length(g.entries)
    n_loci_alleles::Int = length(g.loci_alleles)
    n_loci::Int = 0
    max_n_alleles::Int = 0
    chr::String = ""
    pos::Int = 0
    for locus in g.loci_alleles
        # locus = g.loci_alleles[1]
        locus_ids::Array{String,1} = split(locus, '\t')
        if n_loci == 0
            chr = locus_ids[1]
            pos = parse(Int, locus_ids[2])
            max_n_alleles = length(split(locus_ids[3], '|'))
            n_loci += 1
        else
            if ((chr == locus_ids[1]) && (pos != parse(Int, locus_ids[2]))) ||
               (chr != locus_ids[1])
                chr = locus_ids[1]
                pos = parse(Int, locus_ids[2])
                n_alleles = length(split(locus_ids[3], '|'))
                max_n_alleles = max_n_alleles < n_alleles ? n_alleles : max_n_alleles
                n_loci += 1
            end
        end
    end
    (n_entries, n_loci_alleles, n_loci, max_n_alleles)
end
