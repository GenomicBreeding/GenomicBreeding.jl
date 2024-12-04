"""
Genomes struct containing unique entries and loci where allele frequencies can have missing values

# Fields
- entries: names of the `n` entries or samples
- loci: names of the `p` loci including the chromsome or scaffold name, position, all alleles, and current allele separated by tabs
- allele_frequencies: `n x p` matrix of allele frequencies between 0 and 1 which can have missing values
- mask: `n x p` matrix of boolean mask for selective analyses and slicing

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2, p=2)
Genomes([#undef, #undef], [#undef, #undef], Union{Missing, Real}[#undef #undef; #undef #undef], Bool[0 0; 0 0])

julia> fieldnames(Genomes)
(:entries, :loci, :allele_frequencies, :mask)

julia> genomes.entries = ["entry_1", "entry_2"];

julia> genomes.loci = ["chr1_12345_A|T_A", "chr2_678910_C|D_D"];

julia> genomes.allele_frequencies = [0.5 0.25; 0.9 missing];

julia> genomes.mask = [true true; true false];

julia> genomes
Genomes(["entry_1", "entry_2"], ["chr1_12345_A|T_A", "chr2_678910_C|D_D"], Union{Missing, Real}[0.5 0.25; 0.9 missing], Bool[1 1; 1 0])
```
"""
mutable struct Genomes
    entries::Array{String,1}
    loci::Array{String,1}
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
Check dimension compatibility of the fields of the Genomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> g = Genomes(n=2,p=4);

julia> check(g)
true

julia> g.entries = ["beaking_change"];

julia> check(g)
false
```
"""
function check(g::Genomes)::Bool
    n, p = size(g.allele_frequencies)
    if (n != length(g.entries)) || (p != length(g.loci)) || ((n, p) != size(g.mask))
        return false
    end
    if !isa(g.entries, Array{String,1}) ||
       !isa(g.loci, Array{String,1}) ||
       !isa(g.allele_frequencies, Array{Union{Real,Missing},2}) ||
       !isa(g.mask, Array{Bool,2})
        return false
    end
    return true
end
