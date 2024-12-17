function histallelefreqs(genomes::Genomes)::Nothing
    q::Vector{Float64} = filter(
        !ismissing,
        reshape(genomes.allele_frequencies, (prod(size(genomes.allele_frequencies)), 1)),
    )
    plt_histogram =
        UnicodePlots.histogram(q, vertical = true, xlim = (0.0, 1.0), nbins = 50)
    display(plt_histogram)
    nothing
end


function simulatemating(;
    parent_genomes::Genomes,
    n_generations::Int = 1,
    pop_size_per_gen::Vector{Int64} = [100],
    seed::Int64 = 42,
)::Genomes
    # genomes = simulategenomes(n=1_000, μ_β_params=(0.5,0.5)); parent_genomes = slice(genomes, idx_entries=collect(1:5), idx_loci_alleles=collect(1:length(genomes.loci_alleles))); n_generations::Int = 10; pop_size_per_gen::Vector{Int64} = [100]; seed::Int64 = 42;
    # genomes = simulategenomes(n=1_000, μ_β_params=(2.0,2.0)); parent_genomes = slice(genomes, idx_entries=collect(1:5), idx_loci_alleles=collect(1:length(genomes.loci_alleles))); n_generations::Int = 10; pop_size_per_gen::Vector{Int64} = [100]; seed::Int64 = 42;
    # Check input arguments
    if !checkdims(parent_genomes)
        throw(ArgumentError("Error in the parents' genomes input"))
    end
    if sum(ismissing.(parent_genomes.allele_frequencies)) != 0
        throw(
            ArgumentError(
                "We expect no missing values in the allele frequencies of the parents. Please consider filtering them out or imputing.",
            ),
        )
    end
    # If the vector of population sizes per generation is less than the requested number of generations the we replicate
    if length(pop_size_per_gen) < n_generations
        pop_size_per_gen = repeat(
            pop_size_per_gen,
            outer = Int(ceil(n_generations / length(pop_size_per_gen))),
        )
    end
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Extract loci names
    chromosomes_per_locus_allele::Array{String,1}, _, _ = loci_alleles(parent_genomes)
    chromosomes_per_locus::Array{String,1},
    _,
    loci_ini_idx::Vector{Int64},
    loci_fin_idx::Vector{Int64} = loci(parent_genomes)
    unique_chromosomes::Array{String,1} = unique(chromosomes_per_locus)
    # Iterate across generations
    histallelefreqs(parent_genomes)
    for t = 1:n_generations
        # t = 1
        progeny_genomes =
            Genomes(n = pop_size_per_gen[t], p = length(parent_genomes.loci_alleles))
        progeny_genomes.entries =
            [string("progeny_t", t, "_", i) for i = 1:pop_size_per_gen[t]]
        progeny_genomes.loci_alleles = parent_genomes.loci_alleles
        for chr in unique_chromosomes
            # chr = unique_chromosomes[1]
            idx_loci_alleles::Vector{Int64} = findall(chromosomes_per_locus_allele .== chr)
            idx_loci::Vector{Int64} = findall(chromosomes_per_locus .== chr)
            allele_freqs::Matrix{Float64} =
                parent_genomes.allele_frequencies[:, idx_loci_alleles]
            μ::Vector{Float64} = mean(allele_freqs, dims = 1)[1, :]
            Σ::Matrix{Float64} = StatsBase.cov(allele_freqs)
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
            # Sample the progeny allele frequencies
            progeny_allele_freqs::Matrix{Float64} =
                rand(rng, mvnormal_distribution, pop_size_per_gen[t])'
            # Restrict allele frequencies between zero and one
            progeny_allele_freqs[progeny_allele_freqs.>1.0] .= 1.0
            progeny_allele_freqs[progeny_allele_freqs.<0.0] .= 0.0
            # Make sure allele frequencies sum up to one for multiallelic loci
            for j in eachindex(loci_ini_idx[idx_loci])
                # j = 1
                # println(j)
                idx_ini = loci_ini_idx[idx_loci][j]
                idx_fin = loci_fin_idx[idx_loci][j]
                a = idx_fin - idx_ini
                sum_of_prev_allele_freqs::Vector{Float64} =
                    fill(0.0, size(progeny_allele_freqs, 1))
                for k = j:(j+a)
                    # k = idx_fin
                    sum_of_prev_allele_freqs =
                        sum_of_prev_allele_freqs + progeny_allele_freqs[:, k]
                end
                idx_overloaded = findall(sum_of_prev_allele_freqs .> 1.0)
                # Rescale so that the allele frequncies sum up to one
                progeny_allele_freqs[idx_overloaded, j:(j+a)] ./=
                    sum_of_prev_allele_freqs[idx_overloaded]
            end
            # Update the progenies' genomes
            progeny_genomes.allele_frequencies[:, idx_loci_alleles] = progeny_allele_freqs
        end
        # Update the parent for the next generation
        parent_genomes = progeny_genomes
        println(string("Progenies at generation: ", t))
        histallelefreqs(parent_genomes)
    end
    parent_genomes
end
