"""
# Phenomes struct

Constains unique entries and traits where phenotype data can have missing values

## Constructor

```julia
Phenomes(; n::Int64 = 1, t::Int64 = 2)
```

## Fields
- `entries`: names of the `n` entries or samples
- `populations`: name/s of the population/s each entries or samples belong to
- `traits`: names of the `t` traits
- `phenotypes`: `n x t` matrix of numeric (`R`) phenotype data which can have missing values
- `mask`: `n x t` matrix of boolean mask for selective analyses and slicing

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> phenomes = Phenomes(n=2, t=2)
Phenomes(["", ""], ["", ""], ["", ""], Union{Missing, Float64}[missing missing; missing missing], Bool[0 0; 0 0])

julia> fieldnames(Phenomes)
(:entries, :populations, :traits, :phenotypes, :mask)

julia> phenomes.entries = ["entry_1", "entry_2"];

julia> phenomes.populations = ["pop_A", "pop_B"];

julia> phenomes.traits = ["height", "yield"];

julia> phenomes.phenotypes = [200.0 2.5; 150.0 missing];

julia> phenomes.mask = [true true; true false];

julia> phenomes
Phenomes(["entry_1", "entry_2"], ["pop_A", "pop_B"], ["height", "yield"], Union{Missing, Float64}[200.0 2.5; 150.0 missing], Bool[1 1; 1 0])
```
"""
mutable struct Phenomes
    entries::Array{String,1}
    populations::Array{String,1}
    traits::Array{String,1}
    phenotypes::Array{Union{Float64,Missing},2}
    mask::Matrix{Bool}
    function Phenomes(; n::Int64 = 1, t::Int64 = 2)
        new(fill("", n), fill("", n), fill("", t), fill(missing, n, t), fill(false, n, t))
    end
end

"""
    checkdims(y::Phenomes)::Bool

Check dimension compatibility of the fields of the Phenomes struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> y = Phenomes(n=2, t=2);

julia> checkdims(y)
false

julia> y.entries = ["entry_1", "entry_2"];

julia> y.traits = ["trait_1", "trait_2"];

julia> checkdims(y)
true
```
"""
function checkdims(y::Phenomes)::Bool
    n, p = size(y.phenotypes)
    if (n != length(y.entries)) ||
       (n != length(unique(y.entries))) ||
       (n != length(y.populations)) ||
       (p != length(y.traits)) ||
       (p != length(unique(y.traits))) ||
       ((n, p) != size(y.mask))
        return false
    end
    if !isa(y.entries, Array{String,1}) ||
       !isa(y.populations, Array{String,1}) ||
       !isa(y.traits, Array{String,1}) ||
       !isa(y.phenotypes, Array{Union{Float64,Missing},2}) ||
       !isa(y.mask, Matrix{Bool})
        return false
    end
    return true
end
