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
- `loci_alleles`: names of the `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles) including the 
chromsome or scaffold name, position, all alleles, and current allele separated by tabs ("\\t")
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

julia> genomes.loci_alleles = ["chr1\\t12345\\tA|T\\tA", "chr2\\t678910\\tC|D\\tD"];

julia> genomes.allele_frequencies = [0.5 0.25; 0.9 missing];

julia> genomes.mask = [true true; true false];

julia> genomes
Genomes(["entry_1", "entry_2"], ["pop_1", "pop_1"], ["chr1\\t12345\\tA|T\\tA", "chr2\\t678910\\tC|D\\tD"], Union{Missing, Float64}[0.5 0.25; 0.9 missing], Bool[1 1; 1 0])
```
"""
mutable struct Genomes <: AbstractGB
    entries::Vector{String}
    populations::Vector{String}
    loci_alleles::Vector{String}
    allele_frequencies::Matrix{Union{Float64,Missing}}
    mask::Matrix{Bool}
    function Genomes(; n::Int64 = 1, p::Int64 = 2)
        return new(fill("", n), fill("", n), fill("", p), fill(missing, n, p), fill(false, n, p))
    end
end

"""
    Base.hash(x::Genomes, h::UInt)::UInt

Hash a Genomes struct using the entries, populations and loci_alleles.
We deliberately excluded the allele_frequencies, and mask for efficiency.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = Genomes(n=2, p=2);

julia> typeof(hash(genomes))
UInt64
```
"""
function Base.hash(x::Genomes, h::UInt)::UInt
    # hash(Genomes, hash(x.entries, hash(x.populations, hash(x.loci_alleles, hash(x.allele_frequencies, hash(x.mask, h))))))
    hash(Genomes, hash(x.entries, hash(x.populations, hash(x.loci_alleles, h))))
end


"""
    Base.:(==)(x::Genomes, y::Genomes)::Bool

Equality of Genomes structs using the hash function defined for Genomes structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes_1 = genomes = Genomes(n=2,p=4);

julia> genomes_2 = genomes = Genomes(n=2,p=4);

julia> genomes_3 = genomes = Genomes(n=1,p=2);

julia> genomes_1 == genomes_2
true

julia> genomes_1 == genomes_3
false
```
"""
function Base.:(==)(x::Genomes, y::Genomes)::Bool
    hash(x) == hash(y)
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

julia> genomes.loci_alleles = ["chr1\\t1\\tA|T\\tA", "chr1\\t2\\tC|G\\tG", "chr2\\t3\\tA|T\\tA", "chr2\\t4\\tG|T\\tG"];

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
    if !isa(genomes.entries, Vector{String}) ||
       !isa(genomes.populations, Vector{String}) ||
       !isa(genomes.loci_alleles, Vector{String}) ||
       !isa(genomes.allele_frequencies, Matrix{Union{Float64,Missing}}) ||
       !isa(genomes.mask, Matrix{Bool})
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
        locus_ids::Vector{String} = split(locus, '\t')
        if n_loci == 0
            chr = locus_ids[1]
            pos = parse(Int64, locus_ids[2])
            max_n_alleles = length(split(locus_ids[3], '|'))
            n_chr += 1
            n_loci += 1
        else
            if ((chr == locus_ids[1]) && (pos != parse(Int64, locus_ids[2]))) || (chr != locus_ids[1])
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
    return Dict(
        "n_entries" => n_entries,
        "n_populations" => n_populations,
        "n_loci_alleles" => n_loci_alleles,
        "n_chr" => n_chr,
        "n_loci" => n_loci,
        "max_n_alleles" => max_n_alleles,
    )
end

"""
    loci_alleles(genomes::Genomes)::Tuple{Vector{String},Vector{Int64},Vector{String}}

Extract chromosomes, positions, and alleles across loci-allele combinations

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> chromsomes, positions, alleles = loci_alleles(genomes);

julia> length(chromsomes), length(positions), length(alleles)
(3000, 3000, 3000)
```
"""
function loci_alleles(genomes::Genomes)::Tuple{Vector{String},Vector{Int64},Vector{String}}
    chromosomes::Vector{String} = []
    positions::Vector{Int64} = []
    alleles::Vector{String} = []
    for locus in genomes.loci_alleles
        # locus = genomes.loci_alleles[1]
        locus_ids::Vector{String} = split(locus, '\t')
        push!(chromosomes, locus_ids[1])
        push!(positions, parse(Int64, locus_ids[2]))
        push!(alleles, locus_ids[4])
    end
    return (chromosomes, positions, alleles)
end

"""
    loci(genomes::Genomes)::Tuple{Vector{String},Vector{Int64},Vector{Int64},Vector{Int64}}

Extract chromosome names, positions, start and end indexes of each locus across loci

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> chromsomes, positions, loci_ini_idx, loci_fin_idx = loci(genomes);

julia> length(chromsomes), length(positions), length(loci_ini_idx), length(loci_fin_idx)
(1000, 1000, 1000, 1000)
```
"""
function loci(genomes::Genomes)::Tuple{Vector{String},Vector{Int64},Vector{Int64},Vector{Int64}}
    chromosomes::Vector{String} = []
    positions::Vector{Int64} = []
    loci_ini_idx::Vector{Int64} = []
    loci_fin_idx::Vector{Int64} = []
    idx::Int64 = 0
    for locus in genomes.loci_alleles
        # locus = genomes.loci_alleles[1]
        idx += 1
        locus_ids::Vector{String} = split(locus, '\t')
        if length(chromosomes) == 0
            push!(chromosomes, locus_ids[1])
            push!(positions, parse(Int64, locus_ids[2]))
            push!(loci_ini_idx, idx)
        else
            if ((chromosomes[end] == locus_ids[1]) && (positions[end] != parse(Int64, locus_ids[2]))) ||
               (chromosomes[end] != locus_ids[1])
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
    return (chromosomes, positions, loci_ini_idx, loci_fin_idx)
end

"""
    plot(genomes::Genomes)::Nothing

Plot allele frequencies

# Examples
```
julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

julia> GenomicBreeding.plot(genomes);

```
"""
function plot(genomes::Genomes, seed::Int64 = 42)
    # genomes = simulategenomes(n=100, l=1_000, n_alleles=3, n_populations=3, μ_β_params=(2.0,2.0), verbose=false); seed::Int64=42;
    # Per poulation, using min([250, p]) randomly sampled loci plot:
    #   (1) histogram of allele frequencies per entry,
    #   (2) histogram of mean allele frequencies per locus, and
    #   (3) correlation heatmap of allele frequencies.
    rng::TaskLocalRNG = Random.seed!(seed)
    for pop in unique(genomes.populations)
        # pop = genomes.populations[1]
        p = size(genomes.allele_frequencies, 2)
        idx_row::Vector{Int64} = findall(genomes.populations .== pop)
        idx_col::Vector{Int64} = StatsBase.sample(rng, 1:p, minimum([250, p]); replace = false, ordered = true)
        Q = genomes.allele_frequencies[idx_row, idx_col]
        q::Vector{Float64} = filter(!ismissing, reshape(Q, (length(idx_row) * length(idx_col), 1)))
        plt_1 = UnicodePlots.histogram(
            vcat(q, 1.00 .- q);
            title = string("Per entry allele frequencies (", pop, ")"),
            vertical = true,
            nbins = 50,
        )
        display(plt_1)
        # # Mean allele frequencies across entries unfolded
        μ_q::Vector{Float64} = fill(0.0, p)
        for j = 1:length(idx_col)
            μ_q[j] = mean(skipmissing(Q[:, j]))
        end
        plt_2 = UnicodePlots.histogram(
            vcat(μ_q, 1.00 .- μ_q);
            title = string("Mean allele frequencies (", pop, ")"),
            vertical = true,
            nbins = 50,
        )
        display(plt_2)
        # Correlation between allele frequencies
        idx_col = findall(sum(ismissing.(Q); dims = 1)[1, :] .== 0)
        plt_3 = UnicodePlots.heatmap(
            StatsBase.cor(Q[:, findall(sum(ismissing.(Q); dims = 1)[1, :] .== 0)]);
            height = 50,
            width = 50,
            zlabel = string("Pairwise loci correlation (", pop, ")"),
        )
        display(plt_3)
    end
    return nothing
end

"""
    slice(genomes::Genomes;idx_entries::Vector{Int64},idx_loci_alleles::Vector{Int64})::Genomes

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
function slice(genomes::Genomes; idx_entries::Vector{Int64}, idx_loci_alleles::Vector{Int64})::Genomes
    # genomes::Genomes = simulategenomes(); idx_entries::Vector{Int64}=sample(1:100, 10); idx_loci_alleles::Vector{Int64}=sample(1:10_000, 1000);
    genomes_dims::Dict{String,Int64} = dimensions(genomes)
    n_entries::Int64 = genomes_dims["n_entries"]
    n_loci_alleles::Int64 = genomes_dims["n_loci_alleles"]
    if (minimum(idx_entries) < 1) || (maximum(idx_entries) > n_entries)
        throw(ArgumentError("We accept `idx_entries` from 1 to `n_entries` of `genomes`."))
    end
    if (minimum(idx_loci_alleles) < 1) || (maximum(idx_loci_alleles) > n_loci_alleles)
        throw(ArgumentError("We accept `idx_loci_alleles` from 1 to `n_loci_alleles` of `genomes`."))
    end
    sort!(idx_entries)
    sort!(idx_loci_alleles)
    unique!(idx_entries)
    unique!(idx_loci_alleles)
    n, p = length(idx_entries), length(idx_loci_alleles)
    sliced_genomes::Genomes = Genomes(; n = n, p = p)
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
    return sliced_genomes
end

# """
#     filter(genomes::Genomes; maf::Float64)::Genomes

# Filter Genomes struct by minimum allele frequency

# # Examples
# ```jldoctest; setup = :(using GenomicBreeding)
# julia> genomes = simulategenomes(n=100, l=1_000, n_alleles=4, verbose=false);

# ```
# """
# function filter(genomes::Genomes; maf::Float64)::Genomes
#     # genomes::Genomes = simulategenomes(sparsity=0.01, seed=123456); maf=0.01;
#     if (maf < 0.0) || (maf > 1.0)
#         throw(ArgumentError("We accept `maf` from 0.0 to 1.0."))
#     end
#     q::Array{Float64, 1} = fill(0.0, length(genomes.loci_alleles))
#     for j in eachindex(q)
#         q[j] = mean(skipmissing(genomes.allele_frequencies[:, j]))
#     end
#     idx::Array{Int64, 1} = findall((q .>= maf) .&& (q .<= (1.0-maf)))
#     filtered_genomes::Genomes = slice(genomes, idx_entries=collect(1:length(genomes.entries)); idx_loci_alleles=idx);

#     # Output
#     filtered_genomes
# end
