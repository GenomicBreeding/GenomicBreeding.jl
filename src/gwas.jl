"""
    gwas(input::GBInput)::Vector{String}

Perform genome-wide association study.
Outputs are JLD2 files, each containing the Fit struct for each model-trait-population combination.
Note that the `b_hat` field of each Fit struct contains t-statistics for `gwasols` and z-statistics for the other GWAS models, instead of allele effects in genomic prediction models.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, gwas_models=[gwasols], verbose=false);

julia> fname_test_statistics_jld2s = GenomicBreeding.gwas(input);

julia> length(fname_test_statistics_jld2s) == 3
true
```
"""
function gwas(input::GBInput)::Vector{String}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"], gwas_models=[gwasols])
    # Parse input and prepare the output directory
    populations = input.populations
    traits = input.traits
    models = input.gwas_models
    input.analysis = gwas
    input.fname_out_prefix = prepareoutprefixandoutdir(input)
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = loadgenomesphenomes(input)
    # Instantiate the vector of dataframes and output vector of the resulting filenames where the dataframes will be written into
    if isnothing(populations)
        populations = [nothing] 
    end
    if isnothing(traits)
        traits = phenomes.traits
    end
    model_fits::Vector{Fit} = fill(
        Fit(n = length(genomes.entries), l = length(genomes.loci_alleles)),
        length(populations) * length(traits) * length(models),
    )
    fname_test_statistics_jld2s::Vector{String} = fill("", length(model_fits))
    # Fit the entire data to extract effects per trait per model
    if verbose
        pb = ProgressMeter.Progress(length(model_fits); desc = "Performing GWAS: ")
    end
    for (i, model) in enumerate(models)
        for (j, trait) in enumerate(traits)
            for (k, population) in enumerate(populations) 
                # i = 1; model = models[i]; j = 1; trait = traits[j]; k = 1; population = populations[k]
                idx = ((((i - 1) * length(traits) + j) - 1) * length(populations)) + k
                if verbose
                    println(string("Model: ", model, "| Trait: ", trait, "| Population: ", population))
                end
                # Make sure the output filenames are not problematic by replacing some symbols into underscores
                fname_jld2 =
                    string(fname_out_prefix, "model_", model, "-trait_", trait, "-population_", population, ".jld2")
                problematic_strings::Vector{String} =
                    [" ", "\n", "\t", "(", ")", "&", "|", ":", "=", "+", "*", "%", "@", "!"]
                for s in problematic_strings
                    fname_jld2 = replace(fname_jld2, s => "_")
                end
                if isfile(fname_jld2)
                    throw(
                        ErrorException(
                            "The output file: `" * fname_jld2 * "` exists. Please remove/rename the existing file.",
                        ),
                    )
                end
                Γ, Φ = if isnothing(population)
                    (genomes, slice(phenomes, idx_traits = [j]))
                else
                    idx_entries = findall(genomes.populations .== population)
                    (slice(genomes, idx_entries = idx_entries), slice(phenomes, idx_entries = idx_entries, idx_traits = [j]))
                end
                model_fit = model(genomes = Γ, phenomes = Φ, verbose = false)
                model_fits[idx] = model_fit
                fname_test_statistics_jld2s[idx] = writejld2(model_fits[idx], fname = fname_jld2)
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        ### Correlations between in test-statics
        C = fill(-Inf, length(model_fits), length(model_fits))
        for i in eachindex(model_fits)
            for j in eachindex(model_fits)
                C[i, j] = cor(model_fits[i].b_hat, model_fits[j].b_hat) # unlike GP Fits with an intercept in the first element of b_hat, GWAS Fits do not.
            end
        end
        p = UnicodePlots.heatmap(C, title = "Correlation between allele effects")
        display(p)
    end
    # Output filnames of the tab-delimited files containing the tables of allele effects
    fname_test_statistics_jld2s
end
