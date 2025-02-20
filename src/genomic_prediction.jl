"""
    cv(input::GBInput)::Tuple{Vector{String},Vector{String}}

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

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, fname_out_prefix="GBOutput_cv_3/output-3-", populations=["pop_1", "pop_3"], traits=["trait_1"], n_replications=2, n_folds=3, verbose=false);

julia> fnames_cvs, fnames_notes = cv(input);

julia> length(fnames_cvs) == 32, length(fnames_notes) == 0
(true, true)
```
"""
function cv(input::GBInput)::Tuple{Vector{String},Vector{String}}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input and prepare the output directory
    bulk_cv = input.bulk_cv
    models = input.models
    n_folds = input.n_folds
    n_replications = input.n_replications
    input.analysis = cv
    input.fname_out_prefix = prepareoutprefixandoutdir(input)
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = loadgenomesphenomes(input)
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
    fit(input::GBInput)::Vector{String}

Extract allele effects by fitting the models without cross-validation.
Outputs are JLD2 files, each containing the Fit struct for each model-trait-population combination.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=GenomicBreeding.fit, fname_out_prefix="GBOutput_fit_4/output-4-", populations=["pop_1", "pop_3"], traits=["trait_1"], models=[bayesa, bayesb], verbose=false);

julia> fname_allele_effects_jld2s = GenomicBreeding.fit(input);

julia> length(fname_allele_effects_jld2s) == 4
true
```
"""
function fit(input::GBInput)::Vector{String}
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # Parse input and prepare the output directory
    models = input.models
    n_iter = input.n_iter
    n_burnin = input.n_burnin
    input.analysis = fit
    input.fname_out_prefix = prepareoutprefixandoutdir(input)
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load and merge the genomes and phenomes
    genomes, phenomes = loadgenomesphenomes(input)
    # Define the selection populations and traits
    populations = if !isnothing(input.populations)
        input.populations
    else
        sort(unique(phenomes.populations))
    end
    traits = if !isnothing(input.traits)
        input.traits
    else
        sort(unique(phenomes.traits))
    end
    # Instantiate the vector of dataframes and output vector of the resulting filenames where the dataframes will be written into
    model_fits::Vector{Fit} = []
    fname_allele_effects_jld2s::Vector{String} = []
    # Fit the entire data to extract effects per trait per model
    if verbose
        pb = ProgressMeter.Progress(length(model_fits); desc = "Extracting allele effects: ")
    end
    for model in models
        for trait in traits
            for population in populations
                # model = models[1]; trait = traits[1]; population = populations[1]
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
                Γ, Φ = if population == "bulk"
                    (genomes, phenomes)
                else
                    idx_entries = findall(genomes.populations .== population)
                    (slice(genomes, idx_entries = idx_entries), slice(phenomes, idx_entries = idx_entries))
                end
                model_fit = if !isnothing(match(Regex("bayes", "i"), string(model)))
                    model(genomes = Γ, phenomes = Φ, n_iter = n_iter, n_burnin = n_burnin, verbose = false)
                else
                    model(genomes = Γ, phenomes = Φ, verbose = false)
                end
                push!(model_fits, model_fit)
                push!(fname_allele_effects_jld2s, writejld2(model_fit, fname = fname_jld2))
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        ### Correlations between in allele effects
        C = fill(-Inf, length(model_fits), length(model_fits))
        for i in eachindex(model_fits)
            for j in eachindex(model_fits)
                C[i, j] = cor(model_fits[i].b_hat[2:end], model_fits[j].b_hat[2:end])
                # model1 = model_fits[i].model[1]; trait1 = model_fits[i].trait[1]; population1 = model_fits[i].population[1]
                # model2 = model_fits[j].model[1]; trait2 = model_fits[j].trait[1]; population2 = model_fits[j].population[1]
                # println(join([model1, trait1, population1], "|"), " vs ",join([model2, trait2, population2], "|"), ": ", round(C[i, j], digits=2))
            end
        end
        p = UnicodePlots.heatmap(C, title = "Correlation between allele effects")
        display(p)
    end
    # Output filnames of the tab-delimited files containing the tables of allele effects
    fname_allele_effects_jld2s
end

"""
    predict(input::GBInput)::String

Predict a Phenomes struct, i.e. trait values using the allele effects from `GenomicBreeding.fit(...)` and the input Genomes struct.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, fname_out_prefix="GBOutput_predict_5/output-5-", analysis=GenomicBreeding.fit, verbose=false);

julia> input.fname_allele_effects_jld2s = GenomicBreeding.fit(input);

julia> input.analysis = GenomicBreeding.predict;

julia> fname_phenomes_predicted = GenomicBreeding.predict(input);

julia> phenomes_predicted = readjld2(Phenomes, fname=fname_phenomes_predicted);

julia> dimensions(phenomes_predicted)
Dict{String, Int64} with 8 entries:
  "n_total"       => 5400
  "n_zeroes"      => 0
  "n_nan"         => 0
  "n_entries"     => 300
  "n_traits"      => 18
  "n_inf"         => 0
  "n_populations" => 3
  "n_missing"     => 0
```
"""
function predict(input::GBInput)::String
    # genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # input.fname_allele_effects_jld2s = GenomicBreeding.fit(input)
    # Parse input and prepare the output directory
    input.analysis = predict
    input.fname_out_prefix = prepareoutprefixandoutdir(input)
    fname_out_prefix = input.fname_out_prefix
    fname_allele_effects_jld2s = input.fname_allele_effects_jld2s
    if length(fname_allele_effects_jld2s) == 0
        throw(
            ArgumentError(
                "No model specified for `" *
                "predict" *
                "`. Please specify the Fit struct/s containing the allele effects (i.e. JLD2 output of `GenomicBreeding.fit(...)`) you wish to use.",
            ),
        )
    end
    # Load and merge the genomes and dummy phenomes (as input.fname_pheno is empty, i.e. "")
    genomes, _dummy_phenomes = loadgenomesphenomes(input)
    # Instantiate the predicted Phenomes struct
    n = length(genomes.entries)
    t = length(fname_allele_effects_jld2s)
    phenomes_predicted = Phenomes(n = n, t = t)
    phenomes_predicted.entries = genomes.entries
    phenomes_predicted.populations = genomes.populations
    # Populate the phenotypes matrix and the trait names vector
    model_fits::Vector{Fit} = loadfits(input)
    for (j, (fname_jld2, model_fit)) in enumerate(zip(fname_allele_effects_jld2s, model_fits))
        # j = 1; fname_jld2 = fname_allele_effects_jld2s[j]; model_fit = model_fits[j];
        y_pred = GBModels.predict(fit = model_fit, genomes = genomes, idx_entries = collect(1:n))
        phenomes_predicted.phenotypes[:, j] = y_pred
        phenomes_predicted.traits[j] = replace(basename(fname_jld2), ".jld2" => "")
    end
    # Save the predicted phenomes
    if !checkdims(phenomes_predicted)
        throw(ErrorException("Error computing GEBVs from: " * string(input)))
    end
    fname_jld2 = if length(fname_allele_effects_jld2s) == 1
        string(
            fname_out_prefix,
            replace(basename(fname_allele_effects_jld2s[1]), ".jld2" => ""),
            "-predicted_phenomes.jld2",
        )
    else
        string(fname_out_prefix, hash(input), "-predicted_phenomes.jld2")
    end
    writejld2(phenomes_predicted, fname = fname_jld2)
end
