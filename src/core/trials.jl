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
Trials(Union{Missing, Float64}[missing missing], [#undef, #undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef])

julia> fieldnames(Trials)
(:phenotypes, :traits, :years, :seasons, :harvests, :sites, :replications, :blocks, :rows, :cols, :entries, :populations)
```
"""
mutable struct Trials
    phenotypes::Array{Union{Float64,Missing},2}
    traits::Array{String,1}
    years::Array{String,1}
    seasons::Array{String,1}
    harvests::Array{String,1}
    sites::Array{String,1}
    replications::Array{String,1}
    blocks::Array{String,1}
    rows::Array{String,1}
    cols::Array{String,1}
    entries::Array{String,1}
    populations::Array{String,1}
    function Trials(; n::Int64 = 2, t::Int64 = 2)
        new(
            fill(missing, n, t),
            Array{String,1}(undef, t),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
        )
    end
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
    df_ids::DataFrame = DataFrame(
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
    )
    df_phe::DataFrame = DataFrame(trials.phenotypes, :auto)
    rename!(df_phe, trials.traits)
    df_phe.id = 1:length(trials.years)
    df = innerjoin(df_ids, df_phe, on = :id)
    df
end

function plot(trials::Trials)::Bool
    seed = abs(rand(Int, 1)[1])
    trials, _ = simulatetrials(genomes = simulategenomes(seed = seed), seed = seed)
    df::DataFrame = tabularise(trials)
    cor(trials.phenotypes)
    plt = corrplot(trials.phenotypes)


    false
end
