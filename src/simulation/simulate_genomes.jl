"""
# Simulate genomes

## Arguments
- `n`: number of entries (default = 100)
- `l`: number of loci (default = 10_000)
- `n_chroms`: number of chromosomes (default = 7)
- `n_alleles`: number of alleles per locus (default = 2)
- `max_pos`: total length of the genome in base-pairs (bp) (default = 135_000_000)
- `ld_corr_50perc_kb`: distance in bp at which linkage expressed as correlation between a pair of loci is at 50% (default = 1_000)
- `μ_β_params`: the shape parameters of the Beta distribution from which the mean allele frequencies will be sampled 
  (default = (0.5, 0.5); U-shaped distribution; you may use (2.0, 2.0) for a bell-shaped distribution)
- `sparsity`: Proportion of missing data (default = 0.0)
- `seed`: psuedo-random number generator seed for replicability (default = 42)
- `verbose`: Show progress bar and plot the linkage heatmap into an svg file? (default = true)

## Output
- `Genomes`

## Examples
```jldoctest; setup = :(using GenomicBreeding, StatsBase, Random)
julia> genomes = simulategenomes(n=100, l=10_000, n_alleles=3, verbose=false);

julia> length(genomes.entries)
100

julia> length(genomes.populations)
100

julia> length(genomes.loci_alleles)
20000

julia> size(genomes.allele_frequencies)
(100, 20000)

julia> mean(ismissing.(genomes.allele_frequencies))
0.0

julia> rng::TaskLocalRNG = Random.seed!(123);

julia> idx = StatsBase.sample(rng, range(1, 20_000, step=2), 250, replace = false, ordered = true);

julia> correlations = StatsBase.cor(genomes.allele_frequencies[:, idx]);

julia> correlations[10,10] == 1.00
true

julia> correlations[10,10] > correlations[10,250]
true

julia> genomes = simulategenomes(n=100, l=10_000, n_alleles=3, sparsity=0.25, verbose=false);

julia> mean(ismissing.(genomes.allele_frequencies))
0.25
```
"""
function simulategenomes(;
    n::Int64 = 100,
    n_populations::Int64 = 1,
    l::Int64 = 10_000,
    n_chroms::Int64 = 7,
    n_alleles::Int64 = 2,
    max_pos::Int64 = 135_000_000,
    ld_corr_50perc_kb::Int64 = 1_000,
    μ_β_params::Tuple{Float64,Float64} = (0.5, 0.5),
    sparsity::Float64 = 0.0,
    seed::Int64 = 42,
    verbose::Bool = true,
)::Genomes
    # n::Int64=100; n_populations::Int64 = 2; l::Int64=10_000; n_chroms::Int64=7;n_alleles::Int64=3; max_pos::Int64=135_000_000; ld_corr_50perc_kb::Int64=1_000; seed::Int64=42; μ_β_params::Tuple{Float64, Float64} = (0.5, 0.5); sparsity::Float64 = 0.25; verbose::Bool=true;
    # Parameter checks
    if (n < 1) || (n > 1e9)
        throw(ArgumentError("We accept `n` from 1 to 1 billion."))
    end
    if (n_populations < 1) || (n_populations > n)
        throw(ArgumentError("We accept `n_populations` from 1 to `n`."))
    end
    if (l < 2) || (l > 1e9)
        throw(ArgumentError("We accept `l` from 2 to 1 billion."))
    end
    if (n_chroms < 1) || (n_chroms > 1e6)
        throw(ArgumentError("We accept `n_chroms` from 1 to 1 million."))
    end
    if (n_alleles < 2) || (n_alleles > 5)
        throw(ArgumentError("We accept `n` from 2 to 5, which can be A, T, C, G, and D (for deletion)."))
    end
    if (max_pos < 10) || (max_pos > 160e9)
        throw(ArgumentError("We accept `max_pos` from 10 to 160 billion (genome of *Tmesipteris oblanceolata*)."))
    end
    if (ld_corr_50perc_kb > ceil(max_pos / n_chroms))
        throw(
            ArgumentError(
                "The parameter `ld_corr_50perc_kb` should be less than or equal to `ceil(max_pos/n_chroms)`.",
            ),
        )
    end
    if l > max_pos
        throw(ArgumentError("The parameter `l` should be less than or equal to `max_pos`."))
    end
    if n_chroms > l
        throw(ArgumentError("The parameter `n_chroms` should be less than or equal to `l`."))
    end
    # Instantiate the output struct
    p = l * (n_alleles - 1)
    genomes = Genomes(; n = n, p = p)
    # Instantiate the randomisation
    rng::TaskLocalRNG = Random.seed!(seed)
    # Simulate chromosome lengths and number of loci per chromosome
    l1::Int64 = Int64(floor(max_pos / n_chroms))
    l2::Int64 = Int64(floor(l / n_chroms))
    chrom_lengths::Vector{Int64} = [
        if i < n_chroms
            l1
        elseif l1 * n_chroms < max_pos
            l1 + (max_pos - l1 * n_chroms)
        else
            l1
        end for i = 1:n_chroms
    ]
    chrom_loci_counts::Vector{Int64} = [
        if i < n_chroms
            l2
        elseif l2 * n_chroms < l
            l2 + (l - l2 * n_chroms)
        else
            l2
        end for i = 1:n_chroms
    ]
    # Simulate loci-alleles combinations coordinates
    allele_choices::Vector{String} = ["A", "T", "C", "G", "D"]
    allele_weights::Weights{Float64,Float64,Vector{Float64}} =
        StatsBase.Weights([1.0, 1.0, 1.0, 1.0, 0.1] / sum([1.0, 1.0, 1.0, 1.0, 0.1]))
    positions::Array{Array{Int64},1} = fill(Int64[], n_chroms)
    locus_counter::Int64 = 1
    loci_alleles::Vector{String} = Vector{String}(undef, p)
    for i = 1:n_chroms
        positions[i] = StatsBase.sample(rng, 1:chrom_lengths[i], chrom_loci_counts[i]; replace = false, ordered = true)
        for pos in positions[i]
            all_alleles::Vector{String} =
                StatsBase.sample(rng, allele_choices, allele_weights, n_alleles; replace = false, ordered = false)
            alleles::Vector{String} =
                StatsBase.sample(rng, all_alleles, n_alleles - 1; replace = false, ordered = false)
            for j in eachindex(alleles)
                loci_alleles[locus_counter] = join([string("chrom_", i), pos, join(all_alleles, "|"), alleles[j]], "\t")
                locus_counter += 1
            end
        end
    end
    # Group entries into populations
    population_sizes::Array{Int64} = [Int64(floor(n / n_populations)) for i = 1:n_populations]
    sum(population_sizes) < n ? population_sizes[end] += n - sum(population_sizes) : nothing
    populations::Array{String} = []
    idx_population_groupings::Array{Array{Int64}} = []
    for i = 1:n_populations
        append!(
            populations,
            [string("pop_", lpad(i, length(string(n_populations)), "0")) for _ = 1:population_sizes[i]],
        )
        if i == 1
            push!(idx_population_groupings, collect(1:population_sizes[i]))
        else
            push!(idx_population_groupings, collect((cumsum(population_sizes)[i-1]+1):cumsum(population_sizes)[i]))
        end
    end
    # Simulate allele frequencies with linkage disequillibrium by sampling from a multivariate normal distribution with non-spherical variance-covariance matrix
    allele_frequencies::Matrix{Union{Float64,Missing}} = fill(missing, n, p)
    locus_counter = 1
    if verbose
        pb = ProgressMeter.Progress(n_chroms * n * (n_alleles - 1); desc = "Simulating allele frequencies: ")
    end
    for i = 1:n_chroms
        n_loci::Int64 = chrom_loci_counts[i]
        pos::Vector{Int64} = positions[i]
        Σ::Matrix{Float64} = fill(0.0, (n_loci, n_loci))
        r::Float64 = log(2.0) / ((ld_corr_50perc_kb * 1_000) / chrom_lengths[i]) # from f(x) = 0.5 = 1 / exp(r*x); where x = normalised distance between loci
        for idx1 = 1:n_loci
            for idx2 = 1:n_loci
                dist::Float64 = abs(pos[idx1] - pos[idx2]) / chrom_lengths[i]
                Σ[idx1, idx2] = 1 / exp(r * dist)
            end
        end
        for k = 1:n_populations
            # Sample mean allele frequencies per population
            Beta_distibution = Distributions.Beta(μ_β_params[1], μ_β_params[2])
            μ::Vector{Float64} = rand(rng, Beta_distibution, n_loci)
            # Scale the variance-covariance matrix by the allele frequency means
            # such that the closer to fixation (closer to 0.0 or 1.0) the lower the variance
            idx_greater_than_half::Vector{Int64} = findall(μ .> 0.5)
            σ² = copy(μ)
            σ²[idx_greater_than_half] = 1.00 .- μ[idx_greater_than_half]
            Σ .*= σ² * σ²'
            # Make sure that the variance-covariance matrix is positive definite
            max_iter::Int64 = 10
            iter::Int64 = 1
            while !isposdef(Σ) && (iter < max_iter)
                if iter == 1
                    Σ[diagind(Σ)] .+= 1.0e-12
                end
                Σ[diagind(Σ)] .*= 10.0
                iter += 1
            end
            # Define the multivariate normal distribution
            mvnormal_distribution = Distributions.MvNormal(μ, Σ)
            for a = 1:(n_alleles-1)
                for j in idx_population_groupings[k]
                    idx_ini = ((n_alleles - 1) * (locus_counter - 1)) + a
                    idx_fin = ((n_alleles - 1) * ((locus_counter - 1) + (n_loci - 1))) + a
                    allele_freqs::Vector{Float64} = rand(rng, mvnormal_distribution)
                    if a > 1
                        sum_of_prev_allele_freqs::Vector{Float64} = fill(0.0, n_loci)
                        for ap = 1:(a-1)
                            ap_idx_ini = ((n_alleles - 1) * (locus_counter - 1)) + ap
                            ap_idx_fin = ((n_alleles - 1) * ((locus_counter - 1) + (n_loci - 1))) + ap
                            x = allele_frequencies[j, range(ap_idx_ini, ap_idx_fin; step = (n_alleles - 1))]
                            sum_of_prev_allele_freqs = sum_of_prev_allele_freqs + x
                        end
                        allele_freqs = 1 .- sum_of_prev_allele_freqs
                    end
                    allele_freqs[allele_freqs.>1.0] .= 1.0
                    allele_freqs[allele_freqs.<0.0] .= 0.0
                    allele_frequencies[j, range(idx_ini, idx_fin; step = (n_alleles - 1))] = allele_freqs
                    if verbose
                        ProgressMeter.next!(pb)
                    end
                end # entries per population
            end # alleles
        end # populations
        locus_counter += n_loci
    end
    if verbose
        ProgressMeter.finish!(pb)
    end
    # Populate the output struct
    genomes.entries = ["entry_" * lpad(i, length(string(n)), "0") for i = 1:n]
    genomes.populations = populations
    genomes.loci_alleles = loci_alleles
    genomes.allele_frequencies = allele_frequencies
    genomes.mask = fill(true, (n, p))
    # Simulate sparsity
    if sparsity > 0.0
        idx = sample(rng, 0:((n*p)-1), Int64(round(sparsity * n * p)); replace = false)
        idx_rows = (idx .% n) .+ 1
        idx_cols = Int.(floor.(idx ./ n)) .+ 1
        genomes.allele_frequencies[CartesianIndex.(idx_rows, idx_cols)] .= missing
    end
    if verbose
        # # Plots of population 1 just as an example
        # idx = StatsBase.sample(
        #     rng,
        #     range(1, p, step = (n_alleles - 1)),
        #     250,
        #     replace = false,
        #     ordered = true,
        # )
        # Q = allele_frequencies[1:population_sizes[1], idx]
        # q::Vector{Float64} =
        #     filter(!ismissing, reshape(Q, (population_sizes[1] * length(idx), 1)))
        # append!(q, 1.00 .- q)
        # plt_histogram = UnicodePlots.histogram(
        #     q,
        #     title = "Population 1 allele frequencies",
        #     vertical = true,
        #     xlim = (0.0, 1.0),
        #     nbins = 100,
        # )
        # display(plt_histogram)
        # C = StatsBase.cor(Q[:, findall(sum(ismissing.(Q), dims = 1)[1, :] .== 0)])
        # plt_correlation = UnicodePlots.heatmap(
        #     C,
        #     height = 100,
        #     width = 100,
        #     zlabel = "Pairwise loci correlation",
        # )
        # display(plt_correlation)
        plot(genomes)
    end
    ### Check dimensions
    if !checkdims(genomes)
        throw(DimensionMismatch("Error simulating genomes."))
    end
    # Output
    return (genomes)
end
