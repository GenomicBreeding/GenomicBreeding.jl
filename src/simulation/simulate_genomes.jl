using Random
using StatsBase
using Distributions
using ProgressMeter
using Plots
using Dates
# Plots.backend(:plotly)
Plots.backend(:gr)

function simulategenomes(;
    n::Int = 100,
    p::Int = 10_000,
    n_chroms::Int = 7,
    n_alleles::Int = 2,
    max_pos::Int = 135_000_000,
    ld_corr_50perc_kb::Int = 100_000,
    seed::Int = 42,
    verbose::Bool = true,
)::Genomes
    # n::Int=100;p::Int=10_000;n_chroms::Int=7;n_alleles::Int=2; max_pos::Int=135_000_000; ld_corr_50perc_kb::Int=1e6; seed::Int=42; verbose::Bool=true;
    # Parameter checks
    if (n < 1) || (n > 1e9)
        throw(ArgumentError("We accept `n` between 1 and 1 billion."))
    end
    if (p < 2) || (p > 1e9)
        throw(ArgumentError("We accept `p` between 2 and 1 billion."))
    end
    if (n_chroms < 1) || (n_chroms > 1e6)
        throw(ArgumentError("We accept `n_chroms` between 1 and 1 million."))
    end
    if (n_alleles < 2) || (n_alleles > 5)
        throw(
            ArgumentError(
                "We accept `n` between 2 and 5, which can be A, T, C, G, and D (for deletion).",
            ),
        )
    end
    if (max_pos < 10) || (max_pos > 160e9)
        throw(
            ArgumentError(
                "We accept `max_pos` between 10 and 160 billion (genome of *Tmesipteris oblanceolata*).",
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
    if p > max_pos
        throw(ArgumentError("The parameter `p` should be less than or equal to `max_pos`."))
    end
    if n_chroms > p
        throw(
            ArgumentError("The parameter `n_chroms` should be less than or equal to `p`."),
        )
    end
    # Instantiate the randomisation
    Random.seed!(seed)
    # Simulate chromosome lengths and number of loci per chromosome
    l1::Int = Int(floor(max_pos / n_chroms))
    l2::Int = Int(floor(p / n_chroms))
    chrom_lengths::Array{Int,1} = [
        i < n_chroms ? l1 : l1 * n_chroms < max_pos ? l1 + (max_pos - l1 * n_chroms) : l1 for i = 1:n_chroms
    ]
    chrom_loci_counts::Array{Int,1} = [
        i < n_chroms ? l2 : l2 * n_chroms < p ? l2 + (p - l2 * n_chroms) : l2 for
        i = 1:n_chroms
    ]
    # Simulate loci coordinates
    alleles::Array{String,1} = ["A", "T", "C", "G", "D"]
    allele_weights::Weights{Float64,Float64,Vector{Float64}} =
        StatsBase.Weights([1.0, 1.0, 1.0, 1.0, 0.1] / sum([1.0, 1.0, 1.0, 1.0, 0.1]))
    positions::Array{Array{Int},1} = fill(Int[], n_chroms)
    loci::Array{String,1} = []
    for i = 1:n_chroms
        positions[i] = StatsBase.sample(
            1:chrom_lengths[i],
            chrom_loci_counts[i],
            replace = false,
            ordered = true,
        )
        for pos in positions[i]
            all_alleles::Array{String,1} = StatsBase.sample(
                alleles,
                allele_weights,
                n_alleles,
                replace = false,
                ordered = false,
            )
            allele::String = StatsBase.sample(all_alleles, 1)[1]
            push!(
                loci,
                join([string("chrom_", i), pos, join(all_alleles, "|"), allele], "\t"),
            )
        end
    end
    # Simulate allele frequencies with linkage disequillibrium by sampling from a multivariate normal distribution with non-spherical variance-covariance matrix
    allele_frequencies::Array{Union{Real,Missing},2} = fill(missing, n, p)
    locus_counter::Int = 1
    if verbose
        pb = ProgressMeter.Progress(n_chroms * n, desc = "Simulating allele frequencies: ")
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
        μ::Array{Real,1} = rand(uniform_distibution, n_loci)
        mvnormal_distribution = Distributions.MvNormal(μ, Σ)
        for j = 1:n
            allele_frequencies[j, locus_counter:(locus_counter+n_loci-1)] =
                rand(mvnormal_distribution)
            if verbose
                ProgressMeter.next!(pb)
            end
        end
        locus_counter += n_loci
    end
    if verbose
        ProgressMeter.finish!(pb)
        idx = StatsBase.sample(1:p, 250, replace = false, ordered = true)
        C = StatsBase.cor(allele_frequencies[:, idx])
        plt = Plots.heatmap(C)
        fname_svg =
            "simulated_loci_correlations_heatmap-" *
            Dates.format(now(), "yyyymmddHHMMSS") *
            ".svg"
        Plots.svg(plt, fname_svg)
        println("Please find the correlations heatmap of the simulated loci: " * fname_svg)
    end
    # Output
    genomes = Genomes(
        ["entry_" * lpad(i, length(string(n)), "0") for i = 1:n],
        loci,
        allele_frequencies,
        fill(true, (n, p)),
    )
    if !check(genomes)
        throw(DimensionMismatch("Error simulating genomes."))
    end
    return (genomes)
end


