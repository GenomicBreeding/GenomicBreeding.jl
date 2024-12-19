"""
# Trials struct 
    
Contains phenotype data across years, seasons, harvest, sites, populations, replications, blocks, rows, and columns

## Constructor

```julia
Trials(; n::Int64 = 2, p::Int64 = 2)
```

## Fields
- `phenotypes`: `n x t` matrix of numeric phenotype data which can have missing values
- `traits`: names of the traits `t` traits
- `years`: names of the years corresponding to each row in the phenotype matrix
- `seasons`: names of the seasons corresponding to each row in the phenotype matrix
- `harvests`: names of the harvests corresponding to each row in the phenotype matrix
- `sites`: names of the sites corresponding to each row in the phenotype matrix
- `replications`: names of the replications corresponding to each row in the phenotype matrix
- `blocks`: names of the blocks corresponding to each row in the phenotype matrix
- `rows`: names of the rows corresponding to each row in the phenotype matrix
- `cols`: names of the cols corresponding to each row in the phenotype matrix
- `entries`: names of the entries corresponding to each row in the phenotype matrix
- `populations`: names of the populations corresponding to each row in the phenotype matrix

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2)
Trials(Union{Missing, Float64}[missing missing], ["", ""], [""], [""], [""], [""], [""], [""], [""], [""], [""], [""])

julia> fieldnames(Trials)
(:phenotypes, :traits, :years, :seasons, :harvests, :sites, :replications, :blocks, :rows, :cols, :entries, :populations)
```
"""
mutable struct Trials
    phenotypes::Matrix{Union{Float64,Missing}}
    traits::Vector{String}
    years::Vector{String}
    seasons::Vector{String}
    harvests::Vector{String}
    sites::Vector{String}
    replications::Vector{String}
    blocks::Vector{String}
    rows::Vector{String}
    cols::Vector{String}
    entries::Vector{String}
    populations::Vector{String}
    function Trials(; n::Int64 = 2, t::Int64 = 2)
        return new(
            fill(missing, n, t),
            fill("", t),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
            fill("", n),
        )
    end
end


"""
    Base.hash(x::Trials, h::UInt)::UInt

Hash a Trials struct.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=2, t=2);

julia> hash(trials)
0x0b8a0916b19be2b1
```
"""
function Base.hash(x::Trials, h::UInt)::UInt
    hash(Trials, hash(x.phenotypes, hash(x.traits, hash(x.years, hash(x.seasons, hash(x.harvests, hash(x.sites, hash(x.replications, hash(x.blocks, hash(x.rows, hash(x.cols, hash(x.entries, hash(x.populations, h)))))))))))))
end


"""
    Base.:(==)(x::Trials, y::Trials)::Bool

Equality of Trials structs using the hash function defined for Trials structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials_1 = trials = Trials(n=2, t=4);

julia> trials_2 = trials = Trials(n=2, t=4);

julia> trials_3 = trials = Trials(n=1, t=2);

julia> trials_1 == trials_2
true

julia> trials_1 == trials_3
false
```
"""
function Base.:(==)(x::Trials, y::Trials)::Bool
    hash(x) == hash(y)
end


"""
    checkdims(trials::Trials)::Bool

Check dimension compatibility of the fields of the Trials struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2);

julia> trials.entries = ["entry_1"];

julia> checkdims(trials)
true

julia> trials.entries = ["entering_2_entries", "instead_of_just_1"];

julia> checkdims(trials)
false
```
"""
function checkdims(trials::Trials)::Bool
    n, t = size(trials.phenotypes)
    if (t != length(trials.traits)) ||
       (n != length(trials.years)) ||
       (n != length(trials.seasons)) ||
       (n != length(trials.harvests)) ||
       (n != length(trials.sites)) ||
       (n != length(trials.replications)) ||
       (n != length(trials.blocks)) ||
       (n != length(trials.rows)) ||
       (n != length(trials.cols)) ||
       (n != length(trials.entries)) ||
       (n != length(trials.populations))
        return false
    end
    return true
end

function tabularise(trials::Trials)::DataFrame
    # trials::Trials, _ = simulatetrials(genomes = simulategenomes());
    df_ids::DataFrame = DataFrame(;
        id = 1:length(trials.years),
        years = trials.years,
        seasons = trials.seasons,
        harvests = trials.harvests,
        sites = trials.sites,
        replications = trials.replications,
        blocks = trials.blocks,
        rows = trials.rows,
        cols = trials.cols,
        entries = trials.entries,
        populations = trials.populations,
    )
    df_phe::DataFrame = DataFrame(trials.phenotypes, :auto)
    rename!(df_phe, trials.traits)
    df_phe.id = 1:length(trials.years)
    df = innerjoin(df_ids, df_phe; on = :id)
    return df
end

# # TODO: plot raw data
# function plot(trials::Trials)::Bool
#     # trials, _ = simulatetrials(genomes = simulategenomes())
#     df::DataFrame = tabularise(trials)
#     cor(trials.phenotypes)
#     plt = corrplot(trials.phenotypes)

#     false
# end
