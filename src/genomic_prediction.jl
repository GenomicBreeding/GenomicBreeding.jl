"""
    assess(input::GBInput)::Tuple{Vector{String},Vector{String}}

Assess genomic prediction accuracy via replicated k-fold cross-validation.
Outputs are saved as JLD2 (each containing a CV struct per fold, replication, and trait) and possibly text file/s containing notes describing why some jobs failed.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], n_replications=2, n_folds=3, verbose=false);

julia> fnames_cvs, fnames_notes = assess(input);

julia> length(fnames_cvs) == 32, length(fnames_notes) == 0
(true, true)
```
"""
function assess(input::GBInput)::Tuple{Vector{String},Vector{String}}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input
    bulk_cv = input.bulk_cv
    models = input.models
    n_folds = input.n_folds
    n_replications = input.n_replications
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = load(input)
    dim_genomes = dimensions(genomes)
    # List the cross-validation function/s to use
    cv_functions = if !bulk_cv && (dim_genomes["n_populations"] > 1)
        [cvperpopulation, cvpairwisepopulation, cvleaveonepopulationout]
    else
        [cvbulk]
    end
    # Instantiate the vector of output filenames
    fnames_cvs::Vector{String} = []
    fnames_notes::Vector{String} = []
    # Instantiate the output directory if it does not exist
    if !isdir(dirname(fname_out_prefix))
        try
            mkdir(dirname(fname_out_prefix))
        catch
            throw(ArgumentError("Error creating the output directory: `" * dirname(fname_out_prefix) * "`"))
        end
    end
    # Make sure we do not have any problematic symbols in the output filenames
    directory_name = dirname(fname_out_prefix)
    prefix_name = basename(fname_out_prefix)
    problematic_strings::Vector{String} = [" ", "\n", "\t", "(", ")", "&", "|", ":", "=", "+", "*", "%"]
    for s in problematic_strings
        prefix_name = replace(prefix_name, s => "_")
    end
    if fname_out_prefix != joinpath(directory_name, prefix_name)
        @warn "Modifying the user defined `fname_out_prefix` as it contains problematic symbols, i.e. from `" *
              fname_out_prefix *
              "` into `" *
              joinpath(directory_name, prefix_name) *
              "`."
        fname_out_prefix = joinpath(directory_name, prefix_name)
    end
    # Cross-validate
    for f in cv_functions
        # f = cv_functions[1]
        cvs, notes = f(
            genomes = genomes,
            phenomes = phenomes,
            models = models,
            n_replications = n_replications,
            n_folds = n_folds,
            verbose = verbose,
        )
        for cv in cvs
            # cv = cvs[1];
            fname = string(fname_out_prefix, f, "-", hash(cv), ".jld2")
            writejld2(cv, fname = fname)
            push!(fnames_cvs, fname)
        end
        if length(notes) > 0
            fname = string(fname_out_prefix, f, "-notes.txt")
            open(fname, "w") do file
                for note in notes
                    write(file, note * "\n")
                end
            end
            push!(fnames_notes, fname)
        end
    end
    # Output
    fnames_cvs, fnames_notes
end

"""
    extracteffects(input::GBInput)::Vector{String}

Extract allele effects by fitting the models without cross-validation.
The output are filenames of tab-delimited files containing the allele effects per model, trait and population.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], n_replications=2, n_folds=3, verbose=false);

julia> fname_effects = extracteffects(input);

julia> length(fname_effects) == 6
true
```
"""
function extracteffects(input::GBInput)::Vector{String}
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # Parse input
    models = input.models
    n_iter = input.n_iter
    n_burnin = input.n_burnin
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load and merge the genomes and phenomes
    genomes, phenomes = load(input)
    # Instantiate the vector of dataframes and output vector of the resulting filenames where the dataframes will be written into
    populations = sort(unique(phenomes.populations))
    traits = phenomes.traits
    dfs::Vector{DataFrames.DataFrame} = fill(DataFrame(), (1 + length(populations)) * length(traits) * length(models))
    fnames::Vector{String} = fill("", (1 + length(populations)) * length(traits) * length(models))
    # Instantiate the output directory if it does not exist
    if !isdir(dirname(fname_out_prefix))
        try
            mkdir(dirname(fname_out_prefix))
        catch
            throw(ArgumentError("Error creating the output directory: `" * dirname(fname_out_prefix) * "`"))
        end
    end
    # Make sure we do not have any problematic symbols in the output filenames
    directory_name = dirname(fname_out_prefix)
    prefix_name = basename(fname_out_prefix)
    problematic_strings::Vector{String} = [" ", "\n", "\t", "(", ")", "&", "|", ":", "=", "+", "*", "%"]
    for s in problematic_strings
        prefix_name = replace(prefix_name, s => "_")
    end
    if fname_out_prefix != joinpath(directory_name, prefix_name)
        @warn "Modifying the user defined `fname_out_prefix` as it contains problematic symbols, i.e. from `" *
              fname_out_prefix *
              "` into `" *
              joinpath(directory_name, prefix_name) *
              "`."
        fname_out_prefix = joinpath(directory_name, prefix_name)
    end
    # Fit the entire data to extract effects per trait per model
    if verbose
        pb = ProgressMeter.Progress(length(dfs); desc = "Extracting allele effects: ")
    end
    for (i, model) in enumerate(models)
        for (j, trait) in enumerate(traits)
            for (k, population) in enumerate(vcat("bulk", populations))
                # i = 1; model = models[i]; j = 1; trait = traits[j]; k = 1; population = populations[k]
                idx = ((((i - 1) * length(traits) + j) - 1) * (1 + length(populations))) + k
                if verbose
                    println(string("Model: ", model, "| Trait: ", trait, "| Population: ", population))
                end
                fname = string(fname_out_prefix, "model_", model, "-trait_", trait, "-population_", population, ".tsv")
                for s in problematic_strings
                    fname = replace(fname, s => "_")
                end
                if isfile(fname)
                    throw(ErrorException("The output file: `` exists. Please remove/rename the existing file."))
                end
                Γ, Φ = if population == "bulk"
                    (genomes, phenomes)
                else
                    idx_entries = findall(genomes.populations .== population)
                    (slice(genomes, idx_entries = idx_entries), slice(phenomes, idx_entries = idx_entries))
                end
                fit = if !isnothing(match(Regex("bayes", "i"), string(model)))
                    model(genomes = Γ, phenomes = Φ, n_iter = n_iter, n_burnin = n_burnin, verbose = false)
                else
                    model(genomes = Γ, phenomes = Φ, verbose = false)
                end
                dfs[idx] = DataFrame(
                    model = fit.model,
                    trait = fit.trait,
                    population = population,
                    fit_corr = fit.metrics["cor"],
                    b_hat_labels = fit.b_hat_labels,
                    b_hat = fit.b_hat,
                )
                open(fname, "w") do file
                    write(file, join(names(dfs[idx]), "\t") * "\n")
                    for line in eachrow(dfs[idx])
                        write(file, join(replace.(string.(Vector(line)), "\t" => "-"), "\t") * "\n")
                    end
                end
                fnames[idx] = fname
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        ### Correlations between in allele effects
        C = fill(-Inf, length(dfs), length(dfs))
        for i in eachindex(dfs)
            for j in eachindex(dfs)
                C[i, j] = cor(dfs[i].b_hat[2:end], dfs[j].b_hat[2:end])
                # model1 = dfs[i].model[1]; trait1 = dfs[i].trait[1]; population1 = dfs[i].population[1]
                # model2 = dfs[j].model[1]; trait2 = dfs[j].trait[1]; population2 = dfs[j].population[1]
                # println(join([model1, trait1, population1], "|"), " vs ",join([model2, trait2, population2], "|"), ": ", round(C[i, j], digits=2))
            end
        end
        p = UnicodePlots.heatmap(C, title = "Correlation between allele effects")
        display(p)
    end
    # Output filnames of the tab-delimited files containing the tables of allele effects
    fnames
end
