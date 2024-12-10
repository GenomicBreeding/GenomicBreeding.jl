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
julia> θ::Array{Float64,2} = simulateeffects();

julia> sum(abs.(θ - [-0.0886501800782904; -0.596478483888422])) < 0.00001
true
```
"""
function simulateeffects(;
    p::Int64 = 2,
    q::Int64 = 1,
    λ::Float64 = 1.00,
    seed::Int64 = 42,
)::Array{Float64,2}
    rng::TaskLocalRNG = Random.seed!(seed)
    μ_dist::Exponential = Distributions.Exponential(λ)
    μ::Array{Float64,1} = rand(rng, μ_dist, p)
    Σ::Array{Float64,2} = μ * μ'
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
    X::Array{Float64,2} = rand(rng, dist, q)
    X
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
)::Tuple{Array{Float64,2},Array{Float64,2}}
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
    n::Int64, _n_populations::Int64, p::Int64, l::Int64, max_n_alleles::Int64 =
        dimensions(genomes)
    # Number of loci with additive, dominance, and epistasis allele effects (minimum values of 1, 0, and 0 or 2 (if there is non-zero loci with epistasis then we expect to have at least 2 loci interacting), respectively; Note that these are also imposed in the above arguments checks)
    a::Int64 = Int64(maximum([1, round(l * f_additive)]))
    d::Int64 = Int64(maximum([0, round(l * f_additive * f_dominance)]))
    e::Int64 = Int64(maximum([0, round(l * f_additive * f_epistasis)]))
    if e == 1
        e = 2 # if there is one epistatic locus then it should interact with at least one other locus
    end
    # Instantiate the output vectors
    α::Array{Float64,1} = fill(0.0, p) # additive allele effects of the p loci-allele combinations
    δ::Array{Float64,1} = fill(0.0, p) # dominance allele effects of the p loci-allele combinations
    ξ::Array{Float64,1} = fill(0.0, p) # epistasis allele effects of the p loci-allele combinations
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Define the loci coordinates with non-zero genetic effects
    idx_additive::Array{Int64,1} =
        StatsBase.sample(rng, 1:l, a, replace = false, ordered = true)
    idx_dominance::Array{Int64,1} =
        StatsBase.sample(rng, idx_additive, d, replace = false, ordered = true)
    idx_epistasis::Array{Int64,1} =
        StatsBase.sample(rng, idx_additive, e, replace = false, ordered = true)
    # Sample additive allele effects from a multivariate normal distribution with non-spherical covariance matrix
    # Notes:
    #   - We are simulating effects on max_n_alleles - 1 alleles hence assuming the remaining allele has zero relative effect.
    #   - We are using a begin-end block to modularise the additive allele effects simulation and will do the same for the dominance and epistasis allele effects.
    additive_effects_per_entry::Array{Float64,1} = begin
        # Simulate the additive allele effects
        A::Array{Float64,2} = simulateeffects(p = a, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the additive allele loci
        idx_p_additive::Array{Int64,1} = []
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
    dominance_effects_per_entry::Array{Float64,1} = begin
        # Simulate the dominance allele effects
        D::Array{Float64,2} = simulateeffects(p = d, q = 1, seed = seed)
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
    epistasis_effects_per_entry::Array{Float64,1} = begin
        # Simulate the epistasis allele effects
        E::Array{Float64,2} = simulateeffects(p = e, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the epistasis allele loci
        idx_p_epistasis::Array{Int64,1} = []
        for i = 1:(max_n_alleles-1)
            append!(idx_p_epistasis, (idx_epistasis * (max_n_alleles - 1)) .- (i - 1))
        end
        sort!(idx_p_epistasis)
        # Update the epistasis allele effects
        ξ[idx_p_epistasis] = reshape(E', (e * (max_n_alleles - 1), 1))
        # Simulate the epistasis allele effects as the sum of the products over all possible pairs of epistatic alleles of their allele frequencies, and epistatic sllele effects
        epistasis_per_entry::Array{Float64,1} = fill(0.0, n)
        for i = 1:(length(idx_p_epistasis)-1)
            idx_1 = idx_p_epistasis[i]
            for j = (i+1):length(idx_p_epistasis)
                idx_2 = idx_p_epistasis[j]
                epistasis_per_entry .+=
                    genomes.allele_frequencies[:, idx_1] .*
                    genomes.allele_frequencies[:, idx_2] .* ξ[idx_1] .* ξ[idx_2]
            end
        end
        epistasis_per_entry
    end
    (
        hcat(
            additive_effects_per_entry,
            dominance_effects_per_entry,
            epistasis_effects_per_entry,
        ),
        hcat(α, δ, ξ),
    )
end


"""
# Simulate trials

## Arguments
- `genomes`: Genome struct includes the `n` entries x `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles)
- `f_additive`: proportion of the `l` loci with non-zero additive effects on the phenotype
- `f_dominance`: proportion of the `l*f_additive` additive effects loci with additional dominance effects
- `f_epistasis`: proportion of the `l*f_additive` additive effects loci with additional epistasis effects
- `n_traits`: Number of traits (default = 3)
- `n_years`: Number of years (default = 2)
- `n_seasons`: Number of seasons (default = 4)
- `n_harvests`: Number of harvests (default = 2)
- `n_sites`: Number of sites (default = 4)
- `n_replications`: Number of replications (default = 2)
- `n_blocks`: Number of blocks across the entire field layout (default = missing)
- `n_rows`: Number of rows across the entire field layout (default = missing)
- `n_cols`: Number of columns across the entire field layout (default = missing)
- `proportion_of_variance`: `9` x `n_traits` numeric matrix of scaled/non-scaled proportion of variances allocated to 
   genetic and environmental effects (default = missing; values will be sampled from a uniform distribution
   followed by a biased sample on the first row, i.e. additive effects row).
   The rows correspond to the variance allocated to:
    1. additive genetic effects
    2. dominance genetic effects
    3. epistasis genetic effects
    4. years effects
    5. seasons effects
    6. sites effects
    7. environmental interactions
    8. spatial interactions
    9. GxE interactiions
  
- `seed`: Randomisation seed (default = 42)
- `sparsity`: Proportion of missing data (default = 0.0)
- `verbose`: Show trials simulation progress bar? (default = true)

## Outputs
- `Trials` struct of simulated phenotype data
- Vector of `SimulatedEffects` each corresponding to each trait-year-season-harvest-site-replication combination

## Examples
```jldoctest; setup = :(using GenomicBreeding; using StatsBase)
julia> genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3, verbose=false);

julia> trials::Trials, vector_of_effects::Array{SimulatedEffects,1} = simulatetrials(genomes=genomes, sparsity=0.25, verbose=false);

julia> size(trials.phenotypes)
(12800, 3)

julia> size(trials.traits)
(3,)

julia> unique(trials.entries) == genomes.entries
true

julia> unique(trials.populations) == unique(genomes.populations)
true

julia> abs(mean(ismissing.(trials.phenotypes)) - 0.25) < 0.00001
true
```
"""
function simulatetrials(;
    genomes::Genomes,
    f_additive::Float64 = 0.05,
    f_dominance::Float64 = 0.50,
    f_epistasis::Float64 = 0.25,
    n_traits::Int64 = 3,
    n_years::Int64 = 2,
    n_seasons::Int64 = 4,
    n_harvests::Int64 = 2,
    n_sites::Int64 = 4,
    n_replications::Int64 = 2,
    n_blocks::Union{Missing,Int64} = missing,
    n_rows::Union{Missing,Int64} = missing,
    n_cols::Union{Missing,Int64} = missing,
    proportion_of_variance::Union{Missing,Array{Float64,2}} = missing,
    sparsity::Float64 = 0.0,
    seed::Int64 = 42,
    verbose::Bool = true,
)::Tuple{Trials,Array{SimulatedEffects,1}}
    # genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3, verbose=false); f_additive::Float64 = 0.05; f_dominance::Float64 = 0.25; f_epistasis::Float64 = 0.10; n_traits::Int64 = 3; n_years::Int64 = 2; n_seasons::Int64 = 4; n_harvests::Int64 = 2; n_sites::Int64 = 4; n_replications::Int64 = 2; n_blocks::Union{Missing,Int64} = missing; n_rows::Union{Missing,Int64} = missing; n_cols::Union{Missing,Int64} = missing; proportion_of_variance::Union{Missing,Array{Float64,2}} = missing; sparsity::Float64 = 0.25; seed::Int64 = 42; verbose::Bool = false;
    # Argument checks
    n, field_layout, rng, proportion_of_variance = begin
        if !checkdims(genomes)
            throw(ArgumentError("simulategenomiceffects: error in the genomes input"))
        end
        if (n_traits < 1) || (n_traits > 1e6)
            throw(ArgumentError("We accept `n_traits` from 1 to 1 million."))
        end
        if (n_years < 1) || (n_years > 1e6)
            throw(ArgumentError("We accept `n_years` from 1 to 1 million."))
        end
        if (n_seasons < 1) || (n_seasons > 1e6)
            throw(ArgumentError("We accept `n_seasons` from 1 to 1 million."))
        end
        if (n_harvests < 1) || (n_harvests > 1e6)
            throw(ArgumentError("We accept `n_harvests` from 1 to 1 million."))
        end
        if (n_sites < 1) || (n_sites > 1e6)
            throw(ArgumentError("We accept `n_sites` from 1 to 1 million."))
        end
        if (n_replications < 1) || (n_replications > 1e6)
            throw(ArgumentError("We accept `n_replications` from 1 to 1 million."))
        end
        # Genomes dimensions
        n::Int64, n_populations::Int64, p::Int64, l::Int64, max_n_alleles::Int64 =
            dimensions(genomes)
        # Argument checks (continued...)
        if !ismissing(n_blocks)
            if (n * n_replications % n_blocks) > 0
                throw(
                    ArgumentError(
                        "`n_blocks` does not evenly divide the `n*n_replications`. Please revise or use `missing` for automatic `n_blocks` setting.",
                    ),
                )
            end
        end
        if !ismissing(n_rows)
            if (n * n_replications % n_rows) > 0
                throw(
                    ArgumentError(
                        "`n_rows` does not evenly divide the `n*n_replications`. Please revise or use `missing` for automatic `n_rows` setting.",
                    ),
                )
            end
        end
        if !ismissing(n_cols)
            if (n * n_replications % n_cols) > 0
                throw(
                    ArgumentError(
                        "`n_cols` does not evenly divide the `n*n_replications`. Please revise or use `missing` for automatic `n_cols` setting.",
                    ),
                )
            end
        end
        if !ismissing(n_rows) && !ismissing(n_cols)
            if (n * n_replications) != (n_rows * n_cols)
                throw(
                    ArgumentError(
                        "`n_rows*n_cols` should equal `n*n_replications`. Please revise or use `missing` for automatic `n_rows` and `n_cols` setting.",
                    ),
                )
            end
        end
        # Define the field layout
        # Note that we prefer the number of rows to be less than or equal to the number of columns, and
        # that blocks divide the field along columns
        if ismissing(n_rows) && ismissing(n_cols)
            n_rows = Int64(floor(sqrt(n * n_replications)))
            while ((n * n_replications) % n_rows) > 0
                n_rows -= 1
            end
            n_cols = Int64(n * n_replications / n_rows)
        elseif !ismissing(n_rows) && ismissing(n_cols)
            n_cols = Int64(n * n_replications / n_rows)
        elseif ismissing(n_rows) && !ismissing(n_cols)
            n_rows = Int64(n * n_replications / n_cols)
        end
        if !ismissing(n_blocks)
            if n_cols % n_blocks > 0
                throw(
                    ArgumentError(
                        "`n_blocks` does not evenly divide the `n_cols`. Please revise or use `missing` for automatic `n_blocks` setting.",
                    ),
                )
            end
        else
            n_blocks = n_cols
            while ((n_cols % n_blocks) > 0) || (n_blocks == n_cols)
                n_blocks -= 1
            end
        end
        n_cols_per_block::Int64 = Int64(n_cols / n_blocks)
        field_layout::Array{Int64,2} = fill(0, n_rows * n_cols, 4)
        counter = 1
        counter_entries_per_replication = 1
        replication = 1
        for row = 1:n_rows
            block = 1
            counter_cols_per_block = 1
            for col = 1:n_cols
                field_layout[counter, 1] = replication
                field_layout[counter, 2] = block
                field_layout[counter, 3] = row
                field_layout[counter, 4] = col
                if counter_entries_per_replication == n
                    replication += 1
                    counter_entries_per_replication = 0
                end
                if counter_cols_per_block == n_cols_per_block
                    block += 1
                    counter_cols_per_block = 0
                end
                counter += 1
                counter_entries_per_replication += 1
                counter_cols_per_block += 1
            end
        end
        # Set randomisation seed
        rng::TaskLocalRNG = Random.seed!(seed)
        # Argument checks (continued...)
        if ismissing(proportion_of_variance)
            proportion_of_variance::Union{Missing,Array{Float64,2}} = rand(rng, 9, n_traits)
            proportion_of_variance[1, :] .= sample(rng, 1:10, n_traits)
        end
        if size(proportion_of_variance) != (9, n_traits)
            throw(
                ArgumentError(
                    "We expect 11 x `n_traits` matrix. " *
                    "Each row corresponding to the scaled/non-scaled proportion of variance of " *
                    "\t (1) additive genetic effects\n" *
                    "\t (2) dominance genetic effects\n" *
                    "\t (3) epistasis genetic effects\n" *
                    "\t (4) years effects\n" *
                    "\t (5) seasons effects\n" *
                    "\t (6) sites effects\n" *
                    "\t (7) environmental interactions\n" *
                    "\t (8) spatial interactions\n" *
                    "\t (9) GxE interactions",
                ),
            )
        end
        if sum(proportion_of_variance .>= 0.0) < (9 * n_traits)
            throw(
                ArgumentError(
                    "We accept zero or positive values in `proportion_of_variance` matrix.",
                ),
            )
        end
        if (sparsity < 0.0) || (sparsity > 1.0)
            throw(ArgumentError("We accept `sparsity` from 0.0 to 1.0."))
        end
        (n, field_layout, rng, proportion_of_variance)
    end
    # Instantiate output Trials struct
    n_total::Int64 = n_years * n_seasons * n_harvests * n_sites * n_replications * n
    trials::Trials = Trials(n = n_total, t = n_traits)
    trials.phenotypes = fill(missing, (n_total, n_traits))
    trials.traits = fill("", n_traits)
    trials.years = fill("", n_total)
    trials.seasons = fill("", n_total)
    trials.harvests = fill("", n_total)
    trials.sites = fill("", n_total)
    trials.replications = fill("", n_total)
    trials.blocks = fill("", n_total)
    trials.rows = fill("", n_total)
    trials.cols = fill("", n_total)
    trials.entries = fill("", n_total)
    trials.populations = fill("", n_total)
    # Instantiate the simulated environmental, genetic and various interaction effects
    vector_of_effects::Array{SimulatedEffects,1} = []
    # Scale the `proportion_of_variance` per columns
    proportion_of_variance = proportion_of_variance ./ sum(proportion_of_variance, dims = 1)
    # Simulate effects
    if verbose
        pb = ProgressMeter.Progress(
            n_traits * n_years * n_seasons * n_harvests * n_sites * n_replications,
            desc = "Simulating trial data: ",
        )
    end
    for idx_trait = 1:n_traits
        # Trait name
        trials.traits[idx_trait] = string("trait_", idx_trait)
        # Proportion of variances of the effects
        σ²_additive::Float64,
        σ²_dominance::Float64,
        σ²_epistasis::Float64,
        σ²_year::Float64,
        σ²_season::Float64,
        σ²_site::Float64,
        σ²_environmental_interactions::Float64,
        σ²_spatial_interactions::Float64,
        σ²_GxE_interactions::Float64 = proportion_of_variance[:, idx_trait]
        # Genetic effects, where:
        # - G: genetic effects per entry x additive, dominance, & epistasis
        # - B: allele effects per locus-allele combination x additive, dominance, & epistasis
        G::Array{Float64,2}, B::Array{Float64,2} = simulategenomiceffects(
            genomes = genomes,
            f_additive = f_additive,
            f_dominance = f_dominance,
            f_epistasis = f_epistasis,
            seed = seed,
        )
        # Additive environmental effects
        θ_years::Array{Float64,2} = simulateeffects(p = n_years, seed = seed)
        θ_seasons::Array{Float64,2} = simulateeffects(p = n_seasons, seed = seed)
        θ_sites::Array{Float64,2} = simulateeffects(p = n_sites, seed = seed)
        # Environmental interaction effects
        θ_seasons_x_year::Array{Float64,2} =
            simulateeffects(p = n_seasons, q = n_years, seed = seed)
        θ_harvests_x_season_x_year::Array{Float64,2} =
            simulateeffects(p = n_harvests, q = n_years * n_seasons, seed = seed)
        θ_sites_x_harvest_x_season_x_year::Array{Float64,2} =
            simulateeffects(p = n_sites, q = n_years * n_seasons * n_harvests, seed = seed)
        ### Spatial effects (instantiate here then simulate iteratively below for computational efficiency because the covariance matrices of these are usually big)
        θ_replications_x_site_x_harvest_x_season_x_year::Array{Float64,2} =
            fill(0.0, (n_replications, 1)) # replication nested within year, season, harvest, and site
        θ_blocks_x_site_x_harvest_x_season_x_year::Array{Float64,2} =
            fill(0.0, (n_blocks, 1)) # block nested within year, season, harvest, and site but across replications
        θ_rows_x_site_x_harvest_x_season_x_year::Array{Float64,2} = fill(0.0, (n_rows, 1)) # row nested within year, season, harvest, and site but across replications
        θ_cols_x_site_x_harvest_x_season_x_year::Array{Float64,2} = fill(0.0, (n_cols, 1)) # col nested within year, season, harvest, and site but across replications
        # Allele-by-environment interaction effects to these (instantiate here, then simulate iteratively below)
        n_additive::Int64, n_dominance::Int64, n_epistasis::Int64 = sum(B .!= 0.0, dims = 1)
        idx_additve::Array{Int64,1} = sample(
            rng,
            filter(i -> B[i, 1] != 0.0, 1:size(B, 1)),
            Int64(round(n_additive * σ²_GxE_interactions)),
            replace = false,
        )
        idx_dominance::Array{Int64,1} = sample(
            rng,
            filter(i -> B[i, 2] != 0.0, 1:size(B, 1)),
            Int64(round(n_dominance * σ²_GxE_interactions)),
            replace = false,
        )
        idx_epistasis::Array{Int64,1} = sample(
            rng,
            filter(i -> B[i, 3] != 0.0, 1:size(B, 1)),
            Int64(round(n_epistasis * σ²_GxE_interactions)),
            replace = false,
        )
        θ_additive_allele_x_site_x_harvest_x_season_x_year::Array{Float64,2} =
            fill(0.0, (n, 1))
        θ_dominance_allele_x_site_x_harvest_x_season_x_year::Array{Float64,2} =
            fill(0.0, (n, 1))
        θ_epistasis_allele_x_site_x_harvest_x_season_x_year::Array{Float64,2} =
            fill(0.0, (n, 1))
        # Output row counters
        idx_out_ini::Int64 = 1
        idx_out_fin::Int64 = n
        # Iterate across environments
        for idx_year = 1:n_years
            for idx_season = 1:n_seasons
                for idx_harvest = 1:n_harvests
                    for idx_site = 1:n_sites
                        θ_replications_x_site_x_harvest_x_season_x_year =
                            simulateeffects(p = n_replications, seed = seed)
                        θ_blocks_x_site_x_harvest_x_season_x_year =
                            simulateeffects(p = n_blocks, seed = seed)
                        θ_rows_x_site_x_harvest_x_season_x_year =
                            simulateeffects(p = n_rows, seed = seed)
                        θ_cols_x_site_x_harvest_x_season_x_year =
                            simulateeffects(p = n_cols, seed = seed)
                        θ_additive_allele_x_site_x_harvest_x_season_x_year =
                            genomes.allele_frequencies[:, idx_additve] *
                            simulateeffects(p = length(idx_additve), seed = seed)
                        θ_dominance_allele_x_site_x_harvest_x_season_x_year =
                            genomes.allele_frequencies[:, idx_dominance] *
                            simulateeffects(p = length(idx_dominance), seed = seed)
                        θ_epistasis_allele_x_site_x_harvest_x_season_x_year =
                            genomes.allele_frequencies[:, idx_epistasis] *
                            simulateeffects(p = length(idx_epistasis), seed = seed)
                        for idx_replication = 1:n_replications
                            idx_field_layout::Array{Bool,1} =
                                field_layout[:, 1] .== idx_replication
                            # Randomise the layout of the entries
                            idx_randomised_entries::Array{Int64,1} =
                                StatsBase.sample(rng, 1:n, n, replace = false)
                            # Define the indexes of the complex interaction effects
                            idx_ys::Int = (idx_year - 1) * n_seasons + idx_season
                            idx_ysh::Int = (idx_ys - 1) * n_harvests + idx_harvest
                            # Define the effects
                            effects::SimulatedEffects = SimulatedEffects()
                            effects.id = [
                                trials.traits[idx_trait],
                                string("year_", idx_year),
                                string(
                                    "season_",
                                    lpad(idx_season, length(string(n_seasons)), "0"),
                                ),
                                string(
                                    "harvest_",
                                    lpad(idx_harvest, length(string(n_harvests)), "0"),
                                ),
                                string(
                                    "site_",
                                    lpad(idx_site, length(string(n_sites)), "0"),
                                ),
                                string(
                                    "replication_",
                                    lpad(
                                        idx_replication,
                                        length(string(n_replications)),
                                        "0",
                                    ),
                                ),
                            ]
                            effects.year = θ_years[idx_year] * σ²_year
                            effects.season = θ_seasons[idx_season] * σ²_season
                            effects.site = θ_sites[idx_site] * σ²_site
                            effects.seasons_x_year =
                                θ_seasons_x_year[idx_season, idx_year] *
                                σ²_environmental_interactions
                            effects.harvests_x_season_x_year =
                                θ_harvests_x_season_x_year[idx_harvest, idx_ys] *
                                σ²_environmental_interactions
                            effects.sites_x_harvest_x_season_x_year =
                                θ_sites_x_harvest_x_season_x_year[idx_site, idx_ysh] *
                                σ²_environmental_interactions
                            # Note that the layout of the entries remain constant across traits because we have defined field_layout outside these nested for-loops
                            effects.field_layout = field_layout
                            effects.replications_x_site_x_harvest_x_season_x_year =
                                θ_replications_x_site_x_harvest_x_season_x_year[field_layout[
                                    idx_field_layout,
                                    1,
                                ][idx_randomised_entries]] * σ²_spatial_interactions
                            effects.blocks_x_site_x_harvest_x_season_x_year =
                                θ_blocks_x_site_x_harvest_x_season_x_year[field_layout[
                                    idx_field_layout,
                                    2,
                                ][idx_randomised_entries]] .* σ²_spatial_interactions
                            effects.rows_x_site_x_harvest_x_season_x_year =
                                θ_rows_x_site_x_harvest_x_season_x_year[field_layout[
                                    idx_field_layout,
                                    3,
                                ][idx_randomised_entries]] .* σ²_spatial_interactions
                            effects.cols_x_site_x_harvest_x_season_x_year =
                                θ_cols_x_site_x_harvest_x_season_x_year[field_layout[
                                    idx_field_layout,
                                    4,
                                ][idx_randomised_entries]] .* σ²_spatial_interactions
                            effects.additive_genetic = G[:, 1] * σ²_additive
                            effects.dominance_genetic = G[:, 2] * σ²_dominance
                            effects.epistasis_genetic = G[:, 3] * σ²_epistasis
                            effects.additive_allele_x_site_x_harvest_x_season_x_year =
                                θ_additive_allele_x_site_x_harvest_x_season_x_year[:, 1] .*
                                σ²_GxE_interactions
                            effects.dominance_allele_x_site_x_harvest_x_season_x_year =
                                θ_dominance_allele_x_site_x_harvest_x_season_x_year[:, 1] .*
                                σ²_GxE_interactions
                            effects.epistasis_allele_x_site_x_harvest_x_season_x_year =
                                θ_epistasis_allele_x_site_x_harvest_x_season_x_year[:, 1] .*
                                σ²_GxE_interactions
                            if !checkdims(effects)
                                throw(ErrorException("Error simulating effects."))
                            end
                            # Grow the simulated effects vector
                            push!(vector_of_effects, effects)
                            # Populate the phenotypes matrix
                            trials.phenotypes[idx_out_ini:idx_out_fin, idx_trait] =
                                sum(effects)
                            # Avoid redundancies by assigning the times, sites, and spatial effects IDs once, i.e. on trait 1 only
                            if idx_trait == 1
                                trials.years[idx_out_ini:idx_out_fin] =
                                    repeat([effects.id[1]], outer = n)
                                trials.seasons[idx_out_ini:idx_out_fin] =
                                    repeat([effects.id[2]], outer = n)
                                trials.harvests[idx_out_ini:idx_out_fin] =
                                    repeat([effects.id[3]], outer = n)
                                trials.sites[idx_out_ini:idx_out_fin] =
                                    repeat([effects.id[4]], outer = n)
                                trials.replications[idx_out_ini:idx_out_fin] =
                                    repeat([effects.id[5]], outer = n)
                                trials.blocks[idx_out_ini:idx_out_fin] =
                                    "Block_" .*
                                    lpad.(
                                        field_layout[idx_field_layout, 2][idx_randomised_entries],
                                        length(string(n_blocks)),
                                        "0",
                                    )
                                trials.rows[idx_out_ini:idx_out_fin] =
                                    "Row_" .*
                                    lpad.(
                                        field_layout[idx_field_layout, 3][idx_randomised_entries],
                                        length(string(n_rows)),
                                        "0",
                                    )
                                trials.cols[idx_out_ini:idx_out_fin] =
                                    "Col_" .*
                                    lpad.(
                                        field_layout[idx_field_layout, 4][idx_randomised_entries],
                                        length(string(n_cols)),
                                        "0",
                                    )
                                trials.entries[idx_out_ini:idx_out_fin] = genomes.entries
                                trials.populations[idx_out_ini:idx_out_fin] =
                                    genomes.populations
                            end
                            idx_out_ini += n
                            idx_out_fin = (idx_out_ini - 1) + n
                            if verbose
                                ProgressMeter.next!(pb)
                            end
                        end
                    end
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    # Simulate sparsity
    if sparsity > 0.0
        m::Int64, t::Int64 = size(trials.phenotypes)
        idx = sample(rng, 0:((m*t)-1), Int64(round(sparsity * m * t)), replace = false)
        idx_rows = (idx .% m) .+ 1
        idx_cols = Int.(floor.(idx ./ m)) .+ 1
        trials.phenotypes[CartesianIndex.(idx_rows, idx_cols)] .= missing
    end
    # Output check
    if !checkdims(trials)
        throw(DimensionMismatch("Error simulating genomes."))
    end
    if verbose
        # Show distribution of the traits across the entire simulation
        for j = 1:length(trials.traits)
            # j = 1
            trait = trials.traits[j]
            y = filter(yi -> !ismissing(yi), skipmissing(trials.phenotypes[:, j]))
            println(trait)
            println(
                string(
                    "Sparsity: ",
                    round(100.0 * mean(ismissing.(trials.phenotypes[:, j]))),
                    "%",
                ),
            )
            plt = UnicodePlots.histogram(y)
            display(plt)
        end
    end
    # Output
    (trials, vector_of_effects)
end
