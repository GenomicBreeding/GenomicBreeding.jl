"""
    Base.hash(x::Phenomes, h::UInt)::UInt

Hash a Phenomes struct.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> phenomes = Phenomes(n=2, t=2);

julia> typeof(hash(phenomes))
UInt64
```
"""
function Base.hash(x::Phenomes, h::UInt)::UInt
    hash(Phenomes, hash(x.entries, hash(x.populations, hash(x.traits, hash(x.phenotypes, hash(x.mask, h))))))
end


"""
    Base.:(==)(x::Phenomes, y::Phenomes)::Bool

Equality of Phenomes structs using the hash function defined for Phenomes structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> phenomes_1 = phenomes = Phenomes(n=2, t=4);

julia> phenomes_2 = phenomes = Phenomes(n=2, t=4);

julia> phenomes_3 = phenomes = Phenomes(n=1, t=2);

julia> phenomes_1 == phenomes_2
true

julia> phenomes_1 == phenomes_3
false
```
"""
function Base.:(==)(x::Phenomes, y::Phenomes)::Bool
    hash(x) == hash(y)
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
    if !isa(y.entries, Vector{String}) ||
       !isa(y.populations, Vector{String}) ||
       !isa(y.traits, Vector{String}) ||
       !isa(y.phenotypes, Matrix{Union{Float64,Missing}}) ||
       !isa(y.mask, Matrix{Bool})
        return false
    end
    return true
end


# """
# TODO: filter using mask
# """
# function filter(phenomes::Phenomes)::Phenomes end

# """
# TODO: merge 2 Phenomes structs
# """
# function merge(phenomes::Phenomes, other::Phenomes)::Phenomes end
