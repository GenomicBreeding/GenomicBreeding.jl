"""
Abstract super type for all GenomicBreeding.jl custom types
"""
abstract type AbstractGB end


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
mutable struct Phenomes <: AbstractGB
    entries::Vector{String}
    populations::Vector{String}
    traits::Vector{String}
    phenotypes::Matrix{Union{Float64,Missing}}
    mask::Matrix{Bool}
    function Phenomes(; n::Int64 = 1, t::Int64 = 2)
        return new(fill("", n), fill("", n), fill("", t), fill(missing, n, t), fill(false, n, t))
    end
end


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
mutable struct Trials <: AbstractGB
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
# Trial-estimate breeding values (TEBV) struct 
    
Contains trial-estimated breeding values as generated by `analyse(trials::Trials, ...)`.

## Fields
- `traits`: names of the traits `t` traits
- `formulae`: best-fitting formula for each trait
- `models`: best-fitting linear mixed model for each trait
- `df_BLUEs`: vector of data frames of best linear unbiased estimators or fixed effects table of each best fitting model
- `df_BLUPs`: vector of data frames of best linear unbiased predictors or random effects table of each best fitting model
- `phenomes`: vector of Phenomes structs each containing the breeding values

## Examples
```jldoctest; setup = :(using GenomicBreeding, MixedModels, DataFrames)
julia> tebv = TEBV(traits=[""], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1, t=1)]);

julia> tebv.traits
1-element Vector{String}:
 ""
```
"""
mutable struct TEBV <: AbstractGB
    traits::Vector{String}
    formulae::Vector{String}
    models::Vector{LinearMixedModel{Float64}}
    df_BLUEs::Vector{DataFrame}
    df_BLUPs::Vector{DataFrame}
    phenomes::Vector{Phenomes}
    function TEBV(;
        traits::Vector{String},
        formulae::Vector{String},
        models::Vector{LinearMixedModel{Float64}},
        df_BLUEs::Vector{DataFrame},
        df_BLUPs::Vector{DataFrame},
        phenomes::Vector{Phenomes},
    )
        return new(traits, formulae, models, df_BLUEs, df_BLUPs, phenomes)
    end
end


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
mutable struct SimulatedEffects <: AbstractGB
    id::Vector{String}
    year::Float64
    season::Float64
    site::Float64
    seasons_x_year::Float64
    harvests_x_season_x_year::Float64
    sites_x_harvest_x_season_x_year::Float64
    field_layout::Matrix{Int64}
    replications_x_site_x_harvest_x_season_x_year::Vector{Float64}
    blocks_x_site_x_harvest_x_season_x_year::Vector{Float64}
    rows_x_site_x_harvest_x_season_x_year::Vector{Float64}
    cols_x_site_x_harvest_x_season_x_year::Vector{Float64}
    additive_genetic::Vector{Float64}
    dominance_genetic::Vector{Float64}
    epistasis_genetic::Vector{Float64}
    additive_allele_x_site_x_harvest_x_season_x_year::Vector{Float64}
    dominance_allele_x_site_x_harvest_x_season_x_year::Vector{Float64}
    epistasis_allele_x_site_x_harvest_x_season_x_year::Vector{Float64}
    function SimulatedEffects()
        return new(
            repeat([""]; inner = 6),
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
# Cross-validation struct
"""
mutable struct CV <: AbstractGB
    model::String
    replication::String
    fold::String
    entries::Vector{String}
    populations::Vector{String}
    y_true::Vector{Float64}
    y_pred::Vector{Float64}
end