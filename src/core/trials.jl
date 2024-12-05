"""
# Trials struct 
    
Contains phenotype data across years, seasons, harvest, sites, populations, replications, blocks, rows, and columns

## Constructor

```julia
Trials(; n::Int = 2, p::Int = 2)
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
Trials(Union{Missing, Real}[#undef #undef], [#undef, #undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef])

julia> fieldnames(Trials)
(:phenotypes, :traits, :years, :seasons, :harvests, :sites, :replications, :blocks, :rows, :cols, :entries, :populations)
```
"""
mutable struct Trials
    phenotypes::Array{Union{Real,Missing},2}
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
    function Trials(; n::Int = 2, t::Int = 2)
        new(
            Array{Real,2}(undef, n, t),
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
