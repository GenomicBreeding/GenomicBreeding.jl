using Random
using StatsBase
using Distributions
using ProgressMeter
using Dates
using Plots
using StatsPlots
# Plots.backend(:plotly)
Plots.backend(:gr)

"""
# Simulate genomes

## Arguments
- `n`: number of entries (default = 100)
- `l`: number of loci (default = 10_000)
- `n_chroms`: number of chromosomes (default = 7)
- `n_alleles`: number of alleles per locus (default = 2)
- `max_pos`: total length of the genome in base-pairs (bp) (default = 135_000_000)
- `ld_corr_50perc_kb`: distance in bp at which linkage expressed as correlation between a pair of loci is at 50% (default = 100_000)
- `seed`: psuedo-random number generator seed for replicability (default = 42)
- `verbose`: Show progress bar and plot the linkage heatmap into an svg file? (default = false)

## Output
- `Genomes`

## Examples
```jldoctest; setup = :(using GenomicBreeding, StatsBase, Random)
julia> genomes = simulategenomes(n=100, l=10_000, n_alleles=3);

julia> length(genomes.entries)
100

julia> length(genomes.loci_alleles)
20000

julia> size(genomes.allele_frequencies)
(100, 20000)

julia> rng::TaskLocalRNG = Random.seed!(123);

julia> idx = StatsBase.sample(rng, range(1, 20_000, step=2), 250, replace = false, ordered = true);

julia> correlations = StatsBase.cor(genomes.allele_frequencies[:, idx]);

julia> correlations[10,10] == 1.00
true

julia> correlations[10,10] > correlations[10,250]
true
```
"""
function simulategenomes(;
    n::Int = 100,
    l::Int = 10_000,
    n_chroms::Int = 7,
    n_alleles::Int = 2,
    max_pos::Int = 135_000_000,
    ld_corr_50perc_kb::Int = 100_000,
    seed::Int = 42,
    verbose::Bool = false,
)::Genomes
    # n::Int=100;l::Int=10_000;n_chroms::Int=7;n_alleles::Int=3; max_pos::Int=135_000_000; ld_corr_50perc_kb::Int=1e6; seed::Int=42; verbose::Bool=true;
    # Parameter checks
    if (n < 1) || (n > 1e9)
        throw(ArgumentError("We accept `n` from 1 to 1 billion."))
    end
    if (l < 2) || (l > 1e9)
        throw(ArgumentError("We accept `l` from 2 to 1 billion."))
    end
    if (n_chroms < 1) || (n_chroms > 1e6)
        throw(ArgumentError("We accept `n_chroms` from 1 to 1 million."))
    end
    if (n_alleles < 2) || (n_alleles > 5)
        throw(
            ArgumentError(
                "We accept `n` from 2 to 5, which can be A, T, C, G, and D (for deletion).",
            ),
        )
    end
    if (max_pos < 10) || (max_pos > 160e9)
        throw(
            ArgumentError(
                "We accept `max_pos` from 10 to 160 billion (genome of *Tmesipteris oblanceolata*).",
            ),
        )
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
        throw(
            ArgumentError("The parameter `n_chroms` should be less than or equal to `l`."),
        )
    end
    # Instantiate the output struct
    p = l * (n_alleles - 1)
    genomes = Genomes(n = n, p = p)
    # Instantiate the randomisation
    rng::TaskLocalRNG = Random.seed!(seed)
    # Simulate chromosome lengths and number of loci per chromosome
    l1::Int = Int(floor(max_pos / n_chroms))
    l2::Int = Int(floor(l / n_chroms))
    chrom_lengths::Array{Int,1} = [
        i < n_chroms ? l1 : l1 * n_chroms < max_pos ? l1 + (max_pos - l1 * n_chroms) : l1 for i = 1:n_chroms
    ]
    chrom_loci_counts::Array{Int,1} = [
        i < n_chroms ? l2 : l2 * n_chroms < l ? l2 + (l - l2 * n_chroms) : l2 for
        i = 1:n_chroms
    ]
    # Simulate loci-alleles combinations coordinates
    allele_choices::Array{String,1} = ["A", "T", "C", "G", "D"]
    allele_weights::Weights{Float64,Float64,Vector{Float64}} =
        StatsBase.Weights([1.0, 1.0, 1.0, 1.0, 0.1] / sum([1.0, 1.0, 1.0, 1.0, 0.1]))
    positions::Array{Array{Int},1} = fill(Int[], n_chroms)
    locus_counter::Int = 1
    loci_alleles::Array{String,1} = Array{String,1}(undef, p)
    for i = 1:n_chroms
        positions[i] = StatsBase.sample(
            rng,
            1:chrom_lengths[i],
            chrom_loci_counts[i],
            replace = false,
            ordered = true,
        )
        for pos in positions[i]
            all_alleles::Array{String,1} = StatsBase.sample(
                rng,
                allele_choices,
                allele_weights,
                n_alleles,
                replace = false,
                ordered = false,
            )
            alleles::Array{String,1} = StatsBase.sample(
                rng,
                all_alleles,
                n_alleles - 1,
                replace = false,
                ordered = false,
            )
            for j in eachindex(alleles)
                loci_alleles[locus_counter] = join(
                    [string("chrom_", i), pos, join(all_alleles, "|"), alleles[j]],
                    "\t",
                )
                locus_counter += 1
            end
        end
    end
    # Simulate allele frequencies with linkage disequillibrium by sampling from a multivariate normal distribution with non-spherical variance-covariance matrix
    allele_frequencies::Array{Union{Real,Missing},2} = fill(missing, n, p)
    locus_counter = 1
    if verbose
        pb = ProgressMeter.Progress(
            n_chroms * n * (n_alleles - 1),
            desc = "Simulating allele frequencies: ",
        )
    end
    for i = 1:n_chroms
        n_loci::Int = chrom_loci_counts[i]
        pos::Array{Int,1} = positions[i]
        Σ::Array{Real,2} = fill(0.0, (n_loci, n_loci))
        k::Real = log(2.0) / (ld_corr_50perc_kb / chrom_lengths[i]) # from f(x) = 0.5 = 1 / exp(k*x); where x = normalised distance between loci
        for idx1 = 1:n_loci
            for idx2 = 1:n_loci
                dist::Real = abs(pos[idx1] - pos[idx2]) / chrom_lengths[i]
                Σ[idx1, idx2] = 1 / exp(k * dist)
            end
        end
        uniform_distibution = Distributions.Uniform(0.0, 1.0)
        μ::Array{Real,1} = rand(rng, uniform_distibution, n_loci)
        mvnormal_distribution = Distributions.MvNormal(μ, Σ)
        for a = 1:(n_alleles-1)
            for j = 1:n
                idx_ini = ((n_alleles - 1) * (locus_counter - 1)) + a
                idx_fin = ((n_alleles - 1) * ((locus_counter - 1) + (n_loci - 1))) + a
                allele_freqs::Array{Real,1} = rand(rng, mvnormal_distribution)
                if a > 1
                    sum_of_prev_allele_freqs::Array{Real,1} = fill(0.0, n_loci)
                    for ap = 1:(a-1)
                        ap_idx_ini = ((n_alleles - 1) * (locus_counter - 1)) + ap
                        ap_idx_fin =
                            ((n_alleles - 1) * ((locus_counter - 1) + (n_loci - 1))) + ap
                        x = allele_frequencies[
                            j,
                            range(ap_idx_ini, ap_idx_fin, step = (n_alleles - 1)),
                        ]
                        sum_of_prev_allele_freqs = sum_of_prev_allele_freqs + x
                    end
                    allele_freqs = 1 .- sum_of_prev_allele_freqs
                end
                allele_freqs[allele_freqs.>1.0] .= 1.0
                allele_freqs[allele_freqs.<0.0] .= 0.0
                allele_frequencies[j, range(idx_ini, idx_fin, step = (n_alleles - 1))] =
                    allele_freqs
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
        locus_counter += n_loci
    end
    if verbose
        ProgressMeter.finish!(pb)
        idx = StatsBase.sample(
            rng,
            range(1, p, step = (n_alleles - 1)),
            250,
            replace = false,
            ordered = true,
        )
        C = StatsBase.cor(allele_frequencies[:, idx])
        # plt = StatsPlots.density(skipmissing(reshape(allele_frequencies[:, idx], (n*250, 1))))
        plt = Plots.heatmap(C)
        fname_svg =
            "simulated_loci_correlations_heatmap-" *
            Dates.format(now(), "yyyymmddHHMMSS") *
            ".svg"
        Plots.svg(plt, fname_svg)
        println("Please find the correlations heatmap of the simulated loci: " * fname_svg)
    end
    # Output
    genomes.entries = ["entry_" * lpad(i, length(string(n)), "0") for i = 1:n]
    genomes.loci_alleles = loci_alleles
    genomes.allele_frequencies = allele_frequencies
    genomes.mask = fill(true, (n, p))
    if !checkdims(genomes)
        throw(DimensionMismatch("Error simulating genomes."))
    end
    return (genomes)
end
