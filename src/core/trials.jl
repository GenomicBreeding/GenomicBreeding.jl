"""
Trials struct containing phenotype data across years, seasons, harvest, sites, populations, replications, blocks, rows, and columns

# Fields
- entries: names of the `n` entries or samples
- traits: names of the `t` traits
- phenotypes: `n x t` matrix of numeric (`R`) phenotype data which can have missing values
- mask: `n x t` matrix of boolean mask for selective analyses and slicing

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2)
Trials([#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], Union{Missing, Real}[#undef #undef])

julia> fieldnames(Trials)
(:traits, :years, :seasons, :harvests, :sites, :populations, :replications, :rows, :cols, :entries, :phenotypes)
```
"""
mutable struct Trials
    traits::Array{String,1}
    years::Array{String,1}
    seasons::Array{String,1}
    harvests::Array{String,1}
    sites::Array{String,1}
    populations::Array{String,1}
    replications::Array{String,1}
    rows::Array{String,1}
    cols::Array{String,1}
    entries::Array{String,1}
    phenotypes::Array{Union{Real,Missing},2}
    function Trials(; n::Int = 2, t::Int = 2)
        new(
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
            Array{Real,2}(undef, n, t),
        )
    end
end

"""
Check dimension compatibility of the fields of the Trials struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2);

julia> trials.entries = ["entry_1"];

julia> check(trials)
true

julia> trials.entries = ["entering_2_entries", "instead_of_just_1"];

julia> check(trials)
false
```
"""
function check(trials::Trials)::Bool
    n, t = size(trials.phenotypes)
    if (t != length(trials.traits)) ||
       (n != length(trials.years)) ||
       (n != length(trials.seasons)) ||
       (n != length(trials.harvests)) ||
       (n != length(trials.sites)) ||
       (n != length(trials.populations)) ||
       (n != length(trials.replications)) ||
       (n != length(trials.blocks)) ||
       (n != length(trials.rows)) ||
       (n != length(trials.cols)) ||
       (n != length(trials.entries))
        return false
    end
    return true
end
