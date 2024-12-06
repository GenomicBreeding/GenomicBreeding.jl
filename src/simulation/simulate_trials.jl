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
julia> θ::Array{Real,2} = simulateeffects();

julia> sum(abs.(θ - [-0.0886501800782904; -0.596478483888422])) < 0.00001
true
```
"""
function simulateeffects(; p::Int = 2, q::Int = 1, λ::Real = 1.00, seed::Int = 42)
    rng::TaskLocalRNG = Random.seed!(seed)
    μ_dist::Exponential = Distributions.Exponential(λ)
    μ::Array{Real,1} = rand(rng, μ_dist, p)
    Σ::Array{Real,2} = μ * μ'
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
    X::Array{Real,2} = rand(rng, dist, q)
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
julia> genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3);

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
    f_additive::Real = 0.01,
    f_dominance::Real = 0.10,
    f_epistasis::Real = 0.05,
    seed::Int = 42,
)::Tuple{Array{Real,2},Array{Real,2}}
    # genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3); f_additive::Real = 0.01; f_dominance::Real = 0.10; f_epistasis::Real = 0.05; seed::Int = 42;
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
    n::Int, p::Int, l::Int, max_n_alleles::Int = dimensions(genomes)
    # Number of loci with additive, dominance, and epistasis allele effects (minimum values of 1, 0, and 0 or 2 (if there is non-zero loci with epistasis then we expect to have at least 2 loci interacting), respectively; Note that these are also imposed in the above arguments checks)
    a::Int = Int(maximum([1, round(l * f_additive)]))
    d::Int = Int(maximum([0, round(l * f_additive * f_dominance)]))
    e::Int = Int(maximum([0, round(l * f_additive * f_epistasis)]))
    if e == 1
        e = 2 # if there is one epistatic locus then it should interact with at least one other locus
    end
    # Instantiate the output vectors
    g::Array{Real,1} = fill(0.0, n) # genetic effects of the n entries
    α::Array{Real,1} = fill(0.0, p) # additive allele effects of the p loci-allele combinations
    δ::Array{Real,1} = fill(0.0, p) # dominance allele effects of the p loci-allele combinations
    ξ::Array{Real,1} = fill(0.0, p) # epistasis allele effects of the p loci-allele combinations
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Define the loci coordinates with non-zero genetic effects
    idx_additive::Array{Int,1} =
        StatsBase.sample(rng, 1:l, a, replace = false, ordered = true)
    idx_dominance::Array{Int,1} =
        StatsBase.sample(rng, idx_additive, d, replace = false, ordered = true)
    idx_epistasis::Array{Int,1} =
        StatsBase.sample(rng, idx_additive, e, replace = false, ordered = true)
    # Sample additive allele effects from a multivariate normal distribution with non-spherical covariance matrix
    # Notes:
    #   - We are simulating effects on max_n_alleles - 1 alleles hence assuming the remaining allele has zero relative effect.
    #   - We are using a begin-end block to modularise the additive allele effects simulation and will do the same for the dominance and epistasis allele effects.
    additive_effects_per_entry::Array{Real,1} = begin
        # Simulate the additive allele effects
        A::Array{Real,2} = simulateeffects(p = a, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the additive allele loci
        idx_p_additive = vcat(
            (idx_additive * (max_n_alleles - 1)) .- 1,
            idx_additive * (max_n_alleles - 1),
        )
        sort!(idx_p_additive)
        # Update the additive allele effects
        α[idx_p_additive] = reshape(A', (a * (max_n_alleles - 1), 1))
        # Additive effects per entry
        genomes.allele_frequencies * α
    end
    # Sample dominance allele effects from a multivariate normal distribution with non-spherical covariance matrix
    dominance_effects_per_entry::Array{Real,1} = begin
        # Simulate the dominance allele effects
        D::Array{Real,2} = simulateeffects(p = d, q = 1, seed = seed)
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
    epistasis_effects_per_entry::Array{Real,1} = begin
        # Simulate the epistasis allele effects
        E::Array{Real,2} = simulateeffects(p = e, q = (max_n_alleles - 1), seed = seed)
        # Define the loci-alleles combination indexes corresponding to the epistasis allele loci
        idx_p_epistasis = vcat(
            (idx_epistasis * (max_n_alleles - 1)) .- 1,
            idx_epistasis * (max_n_alleles - 1),
        )
        sort!(idx_p_epistasis)
        # Update the epistasis allele effects
        ξ[idx_p_epistasis] = reshape(E', (e * (max_n_alleles - 1), 1))
        # Simulate the epistasis allele effects as the sum of the products over all possible pairs of epistatic alleles of their allele frequencies, and epistatic sllele effects
        epistasis_per_entry::Array{Real,1} = fill(0.0, n)
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

proportion_of_variance::Array{Real,2} = hcat(
    vcat([0, 0, 0], repeat([1], outer = 9)),
    vcat([9, 0, 0], repeat([1], outer = 9)),
    vcat([1, 1, 1], repeat([0], outer = 9)),
)
# df = DataFrame(hcat(trials.years, trials.seasons, trials.harvests, trials.sites, trials.replications, trials.blocks, trials.rows, trials.cols, trials.entries, trials.populations, trials.phenotypes), :auto)
# column_names::Array{String,1} = [
#     "years",
#     "seasons",
#     "harvests",
#     "sites",
#     "replications",
#     "blocks",
#     "rows",
#     "cols",
#     "entries",
#     "population",
# ]
# append!(column_names, trials.traits)
# rename!(df, column_names)
# df_tmp = groupby(df, column_names[9])
# combine(df_tmp, trials.traits[1] => mean)
# combine(df_tmp, trials.traits[1] => std)
# combine(df_tmp, trials.traits[1] => minimum)
# combine(df_tmp, trials.traits[1] => maximum)

"""
# Simulate trials

## Arguments
TODO

## Output
TODO

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3);

julia> trials::Trials = simulatetrials(genomes=genomes);

julia> size(trials.phenotypes)
(12800, 3)

julia> size(trials.traits)
(3,)
```
"""
function simulatetrials(;
    genomes::Genomes,
    n_traits::Int = 3,
    n_years::Int = 2,
    n_seasons::Int = 4,
    n_harvests::Int = 2,
    n_sites::Int = 4,
    n_replications::Int = 2,
    n_blocks::Union{Missing,Int} = missing,
    n_rows::Union{Missing,Int} = missing,
    n_cols::Union{Missing,Int} = missing,
    proportion_of_variance::Array{Real,2} = proportion_of_variance,
    n_populations::Int = 1,
    seed::Int = 42,
    verbose::Bool = false,
)::Trials
    # genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3); n_traits::Int = 3; n_years::Int = 2; n_seasons::Int = 4; n_harvests::Int = 2; n_sites::Int = 4; n_replications::Int = 2; n_blocks::Union{Missing,Int} = missing; n_rows::Union{Missing,Int} = missing; n_cols::Union{Missing,Int} = missing; proportion_of_variance::Array{Real,2}=hcat(vcat([0,0,0], repeat([1], outer=9)), vcat([9,0,0], repeat([1], outer=9)), vcat([1,1,1], repeat([0], outer=9))); n_populations::Int = 1; seed::Int = 42; verbose::Bool = false;
    # Argument checks
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
    n::Int, p::Int, l::Int, max_n_alleles::Int = dimensions(genomes)
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
        n_rows = Int(floor(sqrt(n * n_replications)))
        while ((n * n_replications) % n_rows) > 0
            n_rows -= 1
        end
        n_cols = Int(n * n_replications / n_rows)
    elseif !ismissing(n_rows) && ismissing(n_cols)
        n_cols = Int(n * n_replications / n_rows)
    elseif ismissing(n_rows) && !ismissing(n_cols)
        n_rows = Int(n * n_replications / n_cols)
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
    n_cols_per_block::Int = Int(n_cols / n_blocks)
    field_layout::Array{Int,2} = fill(0, n_rows * n_cols, 4)
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
    # Argument checks (continued...)
    if size(proportion_of_variance) != (12, n_traits)
        throw(
            ArgumentError(
                "We expect 11 x `n_traits` matrix. " *
                "Each row corresponding to the scaled/non-scaled proportion of variance of " *
                "\t (01) additive genetic effects\n" *
                "\t (02) dominance genetic effects\n" *
                "\t (03) epistasis genetic effects\n" *
                "\t (04) years effects\n" *
                "\t (05) seasons effects\n" *
                "\t (06) harvests effects\n" *
                "\t (07) sites effects\n" *
                "\t (08) replications effects\n" *
                "\t (09) blocks effects\n" *
                "\t (10) rows effects\n" *
                "\t (11) cols effects\n" *
                "\t (12) complex interaction effects",
            ),
        )
    end
    if sum(proportion_of_variance .>= 0.0) < (12 * n_traits)
        throw(
            ArgumentError(
                "We accept zero or positive values in `proportion_of_variance` matrix.",
            ),
        )
    end
    if (n_populations < 1) || (n_populations > n)
        throw(ArgumentError("We accept `n_populations` from 1 to `n`."))
    end
    # Instantiate output Trials struct
    n_total::Int = n_years * n_seasons * n_harvests * n_sites * n_replications * n
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
    # Scale the `proportion_of_variance` per columns
    proportion_of_variance = proportion_of_variance ./ sum(proportion_of_variance, dims = 1)
    # Simulate the genetic effects 
    # G: genetic effects per entry x additive, dominance, & epistasis
    # B: allele effects per locus-allele combination x additive, dominance, & epistasis
    G::Array{Real,2}, B::Array{Real,2} = simulategenomiceffects(
        genomes = genomes,
        f_additive = 0.05,
        f_dominance = 0.75,
        f_epistasis = 0.25,
        seed = seed,
    )
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Simulate population groupings of the entries in the genomes input
    populations::Array{String,1} = fill("", n)
    for i = 1:n
        populations[i] = string(
            "Population_",
            lpad(i, length(string(StatsBase.sample(rng, 1:n_populations, 1))), "0"),
        )
    end
    # Simulate effects
    if verbose
        pb = ProgressMeter.Progress(
            n_traits * n_years * n_seasons * n_harvests * n_sites * n_replications,
            desc = "Simulating trial data: ",
        )
    end
    for idx_trait = 1:n_traits
        # idx_trait = 1
        trials.traits[idx_trait] = string("trait_", idx_trait)
        σ²_additive::Real,
        σ²_dominance::Real,
        σ²_epistasis::Real,
        σ²_year::Real,
        σ²_season::Real,
        σ²_harvest::Real,
        σ²_site::Real,
        σ²_replication::Real,
        σ²_block::Real,
        σ²_row::Real,
        σ²_col::Real,
        σ²_complex_interactions::Real = proportion_of_variance[:, idx_trait]
        # Additive environmental effects
        θ_seasons::Array{Real,2} = fill(0.0, (n_seasons, 1))
        θ_harvests::Array{Real,2} = fill(0.0, (n_harvests, 1))
        θ_sites::Array{Real,2} = fill(0.0, (n_sites, 1))
        θ_replications::Array{Real,2} = fill(0.0, (n_replications, 1))
        θ_blocks::Array{Real,2} = fill(0.0, (n_blocks, 1))
        θ_rows::Array{Real,2} = fill(0.0, (n_rows, 1))
        θ_cols::Array{Real,2} = fill(0.0, (n_cols, 1))
        # Limiting environmental interaction effects to these:
        θ_years::Array{Real,2} = simulateeffects(p = n_years, seed = seed)
        θ_years_x_seasons::Array{Real,2} =
            simulateeffects(p = n_years * n_seasons, seed = seed)
        θ_years_x_seasons_x_harvests::Array{Real,2} =
            simulateeffects(p = n_years * n_seasons * n_harvests, seed = seed)
        θ_years_x_seasons_x_harvests_x_sites::Array{Real,2} =
            simulateeffects(p = n_years * n_seasons * n_harvests * n_sites, seed = seed)
        # Limiting GxE interaction effects to these:
        θ_years_x_entries::Array{Real,2} = simulateeffects(p = n_years * n, seed = seed)
        θ_seasons_x_entries::Array{Real,2} = simulateeffects(p = n_seasons * n, seed = seed)
        θ_harvests_x_entries::Array{Real,2} =
            simulateeffects(p = n_harvests * n, seed = seed)
        θ_sites_x_entries::Array{Real,2} = simulateeffects(p = n_sites * n, seed = seed)
        # Output row counters
        idx_out_ini::Int = 1
        idx_out_fin::Int = n
        # Iterate across environments
        for idx_year = 1:n_years
            θ_seasons = simulateeffects(p = n_seasons, seed = seed)
            for idx_season = 1:n_seasons
                θ_harvests = simulateeffects(p = n_harvests, seed = seed)
                for idx_harvest = 1:n_harvests
                    θ_sites = simulateeffects(p = n_sites, seed = seed)
                    for idx_site = 1:n_sites
                        θ_replications = simulateeffects(p = n_replications, seed = seed)
                        θ_blocks = simulateeffects(p = n_blocks, seed = seed)
                        θ_rows = simulateeffects(p = n_rows, seed = seed)
                        θ_cols = simulateeffects(p = n_cols, seed = seed)
                        for idx_replication = 1:n_replications
                            idx_field_layout::Array{Bool,1} =
                                field_layout[:, 1] .== idx_replication
                            # Randomise the layout of the entries
                            idx_randomised_entries::Array{Int,1} =
                                StatsBase.sample(rng, 1:n, n, replace = false)
                            # Define the indexes of the complex interaction effects
                            idx_ys::Int = (idx_year - 1) * n_seasons + idx_season
                            idx_ysh::Int = (idx_ys - 1) * n_harvests + idx_harvest
                            idx_yshs::Int = (idx_ysh - 1) * n_sites + idx_site
                            # Populate the phenotypes matrix
                            trials.phenotypes[idx_out_ini:idx_out_fin, idx_trait] =
                                G[:, 1] * σ²_additive +
                                G[:, 2] * σ²_dominance +
                                G[:, 3] * σ²_epistasis +
                                θ_replications[field_layout[idx_field_layout, 1][idx_randomised_entries]] *
                                σ²_replication +
                                θ_blocks[field_layout[idx_field_layout, 2][idx_randomised_entries]] *
                                σ²_block +
                                θ_rows[field_layout[idx_field_layout, 3][idx_randomised_entries]] *
                                σ²_row +
                                θ_cols[field_layout[idx_field_layout, 4][idx_randomised_entries]] *
                                σ²_col +
                                (
                                    (
                                        θ_years[idx_year] * σ²_year +
                                        θ_seasons[idx_season] * σ²_season +
                                        θ_harvests[idx_harvest] * σ²_harvest +
                                        θ_sites[idx_site] * σ²_site
                                    ) .+
                                    (
                                        (
                                            θ_years_x_seasons[idx_ys] +
                                            θ_years_x_seasons_x_harvests[idx_ysh] +
                                            θ_years_x_seasons_x_harvests_x_sites[idx_yshs]
                                        ) .+ (
                                            θ_years_x_entries[idx_randomised_entries.+(idx_year-1)*n] +
                                            θ_seasons_x_entries[idx_randomised_entries.+(idx_season-1)*n] +
                                            θ_harvests_x_entries[idx_randomised_entries.+(idx_harvest-1)*n] +
                                            θ_sites_x_entries[idx_randomised_entries.+(idx_site-1)*n]
                                        )
                                    ) * σ²_complex_interactions
                                )
                            trials.years[idx_out_ini:idx_out_fin] =
                                repeat([string("year_", idx_year)], outer = n)
                            trials.seasons[idx_out_ini:idx_out_fin] = repeat(
                                [
                                    string(
                                        "season_",
                                        lpad(idx_season, length(string(n_seasons)), "0"),
                                    ),
                                ],
                                outer = n,
                            )
                            trials.harvests[idx_out_ini:idx_out_fin] = repeat(
                                [
                                    string(
                                        "harvest_",
                                        lpad(idx_harvest, length(string(n_harvests)), "0"),
                                    ),
                                ],
                                outer = n,
                            )
                            trials.sites[idx_out_ini:idx_out_fin] = repeat(
                                [
                                    string(
                                        "site_",
                                        lpad(idx_site, length(string(n_sites)), "0"),
                                    ),
                                ],
                                outer = n,
                            )
                            trials.replications[idx_out_ini:idx_out_fin] = repeat(
                                [
                                    string(
                                        "replication_",
                                        lpad(
                                            idx_replication,
                                            length(string(n_replications)),
                                            "0",
                                        ),
                                    ),
                                ],
                                outer = n,
                            )
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
                            trials.populations[idx_out_ini:idx_out_fin] = populations
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
    # Output check
    if !checkdims(trials)
        throw(DimensionMismatch("Error simulating genomes."))
    end
    trials
end
