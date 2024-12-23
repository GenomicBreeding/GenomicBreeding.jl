"""
    Base.hash(x::SimulatedEffects, h::UInt)::UInt

Hash a SimulatedEffects struct.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> effects = SimulatedEffects();

julia> typeof(hash(effects))
UInt64
```
"""
function Base.hash(x::SimulatedEffects, h::UInt)::UInt
    hash(
        SimulatedEffects,
        hash(
            x.id,
            hash(
                x.year,
                hash(
                    x.season,
                    hash(
                        x.site,
                        hash(
                            x.seasons_x_year,
                            hash(
                                x.harvests_x_season_x_year,
                                hash(
                                    x.sites_x_harvest_x_season_x_year,
                                    hash(
                                        x.field_layout,
                                        hash(
                                            x.replications_x_site_x_harvest_x_season_x_year,
                                            hash(
                                                x.blocks_x_site_x_harvest_x_season_x_year,
                                                hash(
                                                    x.rows_x_site_x_harvest_x_season_x_year,
                                                    hash(
                                                        x.cols_x_site_x_harvest_x_season_x_year,
                                                        hash(
                                                            x.additive_genetic,
                                                            hash(
                                                                x.dominance_genetic,
                                                                hash(
                                                                    x.epistasis_genetic,
                                                                    hash(
                                                                        x.additive_allele_x_site_x_harvest_x_season_x_year,
                                                                        hash(
                                                                            x.dominance_allele_x_site_x_harvest_x_season_x_year,
                                                                            hash(
                                                                                x.epistasis_allele_x_site_x_harvest_x_season_x_year,
                                                                                h,
                                                                            ),
                                                                        ),
                                                                    ),
                                                                ),
                                                            ),
                                                        ),
                                                    ),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
end


"""
    Base.:(==)(x::SimulatedEffects, y::SimulatedEffects)::Bool

Equality of SimulatedEffects structs using the hash function defined for SimulatedEffects structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> effects_1 = SimulatedEffects();

julia> effects_2 = SimulatedEffects();

julia> effects_3 = SimulatedEffects(); effects_3.id[1] = "SOMETHING_ELSE";

julia> effects_1 == effects_2
true

julia> effects_1 == effects_3
false
```
"""
function Base.:(==)(x::SimulatedEffects, y::SimulatedEffects)::Bool
    hash(x) == hash(y)
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
function Base.sum(effects::SimulatedEffects)::Vector{Float64}
    ϕ::Vector{Float64} = fill(0.0, size(effects.additive_genetic))
    for name in fieldnames(SimulatedEffects)
        if (name == :id) || (name == :field_layout)
            continue
        end
        ϕ .+= getproperty(effects, name)
    end
    return ϕ
end

"""
# Simulate effects

Sample `p` x `q` effects from a multivariate normal distribution with
`μ~Exp(λ)` and `Σ=μμ'`

## Arguments
- `p`: number of correlated effects to simulate (default = 2)
- `q`: number times to simulate the correlated effects from the same distribution (default = 1)
- `λ`: parameter of the exponential distritbution from which the means will be sampled from (default = 1.00)
- `seed`: randomisation seed (default = 42)

## Output
- `p` x `q` matrix of correlated effects

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> θ::Matrix{Float64} = simulateeffects();

julia> sum(abs.(θ - [-0.0886501800782904; -0.596478483888422])) < 0.00001
true
```
"""
function simulateeffects(; p::Int64 = 2, q::Int64 = 1, λ::Float64 = 1.00, seed::Int64 = 42)::Matrix{Float64}
    # p::Int64 = 20; q::Int64 = 1; λ::Float64 = 1.00; seed::Int64 = 42;
    rng::TaskLocalRNG = Random.seed!(seed)
    μ_dist::Exponential = Distributions.Exponential(λ)
    μ::Vector{Float64} = rand(rng, μ_dist, p)
    # idx_negative_μ::Vector{Int64} = sample(rng, 1:p, Int(floor(rand(rng)*p)))
    # μ[idx_negative_μ] .*= -1.00
    Σ::Matrix{Float64} = μ * μ'
    while abs(det(Σ)) < 1e-12
        Σ[diagind(Σ)] .+= 1.0
        if abs(det(Σ)) >= 1e-12
            break
        else
            μ = rand(rng, μ_dist, p)
            Σ = μ * μ'
        end
    end
    dist::MvNormal = Distributions.MvNormal(μ, Σ)
    X::Matrix{Float64} = rand(rng, dist, q)
    return X
end

"""
# Simulate genomic effects

Simulate additive, dominance, and epistatic effects

## Arguments
- `genomes`: Genome struct includes the `n` entries x `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles)
- `f_additive`: proportion of the `l` loci with non-zero additive effects on the phenotype
- `f_dominance`: proportion of the `l*f_additive` additive effects loci with additional dominance effects
- `f_epistasis`: proportion of the `l*f_additive` additive effects loci with additional epistasis effects

## Outputs
- `n` x `3` matrix of additive, dominance and epistasis effects per entry
- `p` x `3` matrix of additive, dominance and epistasis effects per locus-allele combination

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3, verbose=false);

julia> G, B = simulategenomiceffects(genomes=genomes, f_additive=0.05, f_dominance=0.75, f_epistasis=0.25);

julia> size.([G, B])
2-element Vector{Tuple{Int64, Int64}}:
 (100, 3)
 (4000, 3)

julia> sum(B .!= 0.0, dims=1)
1×3 Matrix{Int64}:
 200  75  50
```

## Details
The additive, dominance, and epistasis allele effects share a common exponential distribution (`λ=1`) from which 
the mean of the effects (`μ`) are sampled, and the covariance matrix is derived (`Σ = μ * μ'`; 
where if `det(Σ)≈0` then we iteratively add 1.00 to the diagonals until it becomes invertible or 10 iterations 
finishes and throws an error). The non-additive or epistasis allele effects were simulated by multiplying the allele 
frequencies of all possible unique pairs of epistasis alleles and their effects.
"""
function simulategenomiceffects(;
    genomes::Genomes,
    f_additive::Float64 = 0.01,
    f_dominance::Float64 = 0.10,
    f_epistasis::Float64 = 0.05,
    seed::Int64 = 42,
)::Tuple{Matrix{Float64},Matrix{Float64}}
    # genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3, verbose=false); f_additive::Float64 = 0.01; f_dominance::Float64 = 0.10; f_epistasis::Float64 = 0.05; seed::Int64 = 42;
    # Argument checks
    if !checkdims(genomes)
        throw(ArgumentError("simulategenomiceffects: error in the genomes input"))
    end
    if (f_additive < 0.0) || (f_additive > 1.0)
        throw(ArgumentError("We accept `f_additive` from 0.00 to 1.00."))
    end
    if (f_dominance < 0.0) || (f_dominance > 1.0)
        throw(ArgumentError("We accept `f_dominance` from 0.00 to 1.00."))
    end
    if (f_epistasis < 0.0) || (f_epistasis > 1.0)
        throw(ArgumentError("We accept `f_epistasis` from 0.00 to 1.00."))
    end
    # Genomes dimensions
    genomes_dims::Dict{String,Int64} = dimensions(genomes)
    n::Int64 = genomes_dims["n_entries"]
    n_populations::Int64 = genomes_dims["n_populations"]
    p::Int64 = genomes_dims["n_loci_alleles"]
    l::Int64 = genomes_dims["n_loci"]
    max_n_alleles::Int64 = genomes_dims["max_n_alleles"]
    # Number of loci with additive, dominance, and epistasis allele effects (minimum values of 1, 0, and 0 or 2 (if there is non-zero loci with epistasis then we expect to have at least 2 loci interacting), respectively; Note that these are also imposed in the above arguments checks)
    a::Int64 = Int64(maximum([1, round(l * f_additive)]))
    d::Int64 = Int64(maximum([0, round(l * f_additive * f_dominance)]))
    e::Int64 = Int64(maximum([0, round(l * f_additive * f_epistasis)]))
    if e == 1
        e = 2 # if there is one epistatic locus then it should interact with at least one other locus
    end
    # Instantiate the output vectors
    α::Vector{Float64} = fill(0.0, p) # additive allele effects of the p loci-allele combinations
    δ::Vector{Float64} = fill(0.0, p) # dominance allele effects of the p loci-allele combinations
    ξ::Vector{Float64} = fill(0.0, p) # epistasis allele effects of the p loci-allele combinations
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Define the loci coordinates with non-zero genetic effects
    idx_additive::Vector{Int64} = StatsBase.sample(rng, 1:l, a; replace = false, ordered = true)
    idx_dominance::Vector{Int64} = StatsBase.sample(rng, idx_additive, d; replace = false, ordered = true)
    idx_epistasis::Vector{Int64} = StatsBase.sample(rng, idx_additive, e; replace = false, ordered = true)
    # Sample additive allele effects from a multivariate normal distribution with non-spherical covariance matrix
    # Notes:
    #   - We are simulating effects on max_n_alleles - 1 alleles hence assuming the remaining allele has zero relative effect.
    #   - We are using a begin-end block to modularise the additive allele effects simulation and will do the same for the dominance and epistasis allele effects.
    additive_effects_per_entry::Vector{Float64} = begin
        # Simulate the additive allele effects
        A::Matrix{Float64} = simulateeffects(; p = a, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the additive allele loci
        idx_p_additive::Vector{Int64} = []
        for i = 1:(max_n_alleles-1)
            append!(idx_p_additive, (idx_additive * (max_n_alleles - 1)) .- (i - 1))
        end
        sort!(idx_p_additive)
        # Update the additive allele effects
        α[idx_p_additive] = reshape(A', (a * (max_n_alleles - 1), 1))
        # Additive effects per entry
        genomes.allele_frequencies * α
    end
    # Sample dominance allele effects from a multivariate normal distribution with non-spherical covariance matrix
    dominance_effects_per_entry::Vector{Float64} = begin
        # Simulate the dominance allele effects
        D::Matrix{Float64} = simulateeffects(; p = d, q = 1, seed = seed)
        # Define the loci-alleles combination indexes corresponding to the first allele per locus with a dominance effect
        idx_p_dominance = (idx_dominance * (max_n_alleles - 1)) .- 1
        sort!(idx_p_dominance)
        # Update the dominance allele effects
        δ[idx_p_dominance] = D[:, 1]
        # Dominance effects per entry
        genomes.allele_frequencies * δ
    end
    # Sample epistasis allele effects from a multivariate normal distribution with non-spherical covariance matrix
    # Notes:
    #   - We are simulating effects on max_n_alleles - 1 alleles hence assuming the remaining allele has zero relative effect.
    #   - Then we simulate the non-additive or epistasis allele effects by multiplying the allele frequencies of 2 epistasis loci and their effects.
    epistasis_effects_per_entry::Vector{Float64} = begin
        # Simulate the epistasis allele effects
        E::Matrix{Float64} = simulateeffects(; p = e, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the epistasis allele loci
        idx_p_epistasis::Vector{Int64} = []
        for i = 1:(max_n_alleles-1)
            append!(idx_p_epistasis, (idx_epistasis * (max_n_alleles - 1)) .- (i - 1))
        end
        sort!(idx_p_epistasis)
        # Update the epistasis allele effects
        ξ[idx_p_epistasis] = reshape(E', (e * (max_n_alleles - 1), 1))
        # Simulate the epistasis allele effects as the sum of the products over all possible pairs of epistatic alleles of their allele frequencies, and epistatic sllele effects
        epistasis_per_entry::Vector{Float64} = fill(0.0, n)
        for i = 1:(length(idx_p_epistasis)-1)
            idx_1 = idx_p_epistasis[i]
            for j = (i+1):length(idx_p_epistasis)
                idx_2 = idx_p_epistasis[j]
                epistasis_per_entry .+=
                    genomes.allele_frequencies[:, idx_1] .* genomes.allele_frequencies[:, idx_2] .* ξ[idx_1] .* ξ[idx_2]
            end
        end
        epistasis_per_entry
    end
    return (hcat(additive_effects_per_entry, dominance_effects_per_entry, epistasis_effects_per_entry), hcat(α, δ, ξ))
end
