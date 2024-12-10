"""
# SimulatedEffects struct 

Contains:
1. Identification, i.e. trait, year, season, harvest, site, and replication
2. Additive environmental effects:
    - year
    - season
    - site
3. Environmental interaction effects:
    - seasons_x_year
    - harvests_x_season_x_year
    - sites_x_harvest_x_season_x_year
4. Spatial effects including the field layout per year-season-harvest-site combination
    - field_layout
    - replications_x_site_x_harvest_x_season_x_year
    - blocks_x_site_x_harvest_x_season_x_year
    - rows_x_site_x_harvest_x_season_x_year
    - cols_x_site_x_harvest_x_season_x_year
5. Genetic effects
    - additive_genetic
    - dominance_genetic
    - epistasis_genetic
6. GxE effects
    - additive_allele_x_site_x_harvest_x_season_x_year
    - dominance_allele_x_site_x_harvest_x_season_x_year
    - epistasis_allele_x_site_x_harvest_x_season_x_year

## Constructor

```julia
SimulatedEffects()
```
"""
mutable struct SimulatedEffects
    id::Array{String,1}
    year::Float64
    season::Float64
    site::Float64
    seasons_x_year::Float64
    harvests_x_season_x_year::Float64
    sites_x_harvest_x_season_x_year::Float64
    field_layout::Array{Int64,2}
    replications_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    blocks_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    rows_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    cols_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    additive_genetic::Array{Float64,1}
    dominance_genetic::Array{Float64,1}
    epistasis_genetic::Array{Float64,1}
    additive_allele_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    dominance_allele_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    epistasis_allele_x_site_x_harvest_x_season_x_year::Array{Float64,1}
    function SimulatedEffects()
        new(
            repeat([""], inner = 6),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            [0 0 0 0; 0 0 0 0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
        )
    end
end


"""
    checkdims(effects::SimulatedEffects)::Bool

Check dimension compatibility of the fields of the SimulatedEffects struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> effects = SimulatedEffects();

julia> checkdims(effects)
true

julia> effects.id = ["beaking_change"];

julia> checkdims(effects)
false
```
"""
function checkdims(effects::SimulatedEffects)::Bool
    n::Int64 = length(effects.replications_x_site_x_harvest_x_season_x_year)
    if (length(effects.id) != 6) ||
       (size(effects.field_layout, 2) != 4) ||
       (n != length(effects.blocks_x_site_x_harvest_x_season_x_year)) ||
       (n != length(effects.rows_x_site_x_harvest_x_season_x_year)) ||
       (n != length(effects.cols_x_site_x_harvest_x_season_x_year)) ||
       (n != length(effects.additive_genetic)) ||
       (n != length(effects.dominance_genetic)) ||
       (n != length(effects.epistasis_genetic)) ||
       (n != length(effects.additive_allele_x_site_x_harvest_x_season_x_year)) ||
       (n != length(effects.dominance_allele_x_site_x_harvest_x_season_x_year)) ||
       (n != length(effects.epistasis_allele_x_site_x_harvest_x_season_x_year))
        return false
    end
    return true
end


"""
    sum(effects::SimulatedEffects)::Tuple{Int64, Int64, Int64}

Sum up the simulated effects to generate the simulated phenotype values

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> effects = SimulatedEffects();

julia> sum(effects)
1-element Vector{Float64}:
 0.0

julia> effects.additive_genetic[1] = pi;

julia> sum(effects)
1-element Vector{Float64}:
 3.141592653589793
```
"""
function Base.sum(effects::SimulatedEffects)::Array{Float64,1}
    ϕ::Array{Float64,1} = fill(0.0, size(effects.additive_genetic))
    for name in fieldnames(SimulatedEffects)
        if (name == :id) || (name == :field_layout)
            continue
        end
        ϕ .+= getproperty(effects, name)
    end
    ϕ
end
