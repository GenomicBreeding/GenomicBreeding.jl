"""
    cv(input::GBInput)::Tuple{Vector{String},Vector{String}}

Assess genomic prediction accuracy via replicated k-fold cross-validation.

# Arguments
- `input::GBInput`: A GBInput struct containing configuration parameters including:
  - `bulk_cv`: Boolean flag for bulk cross-validation
  - `populations`: Vector of population names to analyze
  - `models`: Statistical models to use for prediction
  - `n_folds`: Number of folds for cross-validation
  - `n_replications`: Number of replications for cross-validation
  - `fname_out_prefix`: Prefix for output filenames
  - `verbose`: Boolean flag for detailed output

# Returns
- `Tuple{Vector{String},Vector{String}}`: A tuple containing:
  - First element: Vector of paths to JLD2 files containing CV results
  - Second element: Vector of paths to text files containing error notes

# Details
The function supports three types of cross-validation:
- Single population CV when one population is specified
- Pairwise population CV when two populations are specified
- Both pairwise and leave-one-population-out CV when more than two populations are specified

Results are saved as JLD2 files (one per fold, replication, and trait) and optional text files 
containing notes about failed jobs.

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, bulk_cv=false, fname_out_prefix="GBOutput_cv_3/output-3-", populations=["pop_1", "pop_3"], traits=["trait_1"], n_replications=2, n_folds=3, verbose=false);

julia> fnames_cvs, fnames_notes = cv(input);

julia> length(fnames_cvs) == 4, length(fnames_notes) == 0
(true, true)
```
"""
function cv(input::GBInput)::Tuple{Vector{String},Vector{String}}
    # genomes = GenomicBreedingCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input and prepare the output directory
    bulk_cv = input.bulk_cv
    populations = input.populations
    models = input.models
    n_folds = input.n_folds
    n_replications = input.n_replications
    input.analysis = cv
    input.fname_out_prefix = prepareoutprefixandoutdir(input)
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = loadgenomesphenomes(input)
    # List the cross-validation function/s to use
    cv_functions = if !bulk_cv && !isnothing(populations)
        if length(populations) == 1
            [cvperpopulation]
        elseif length(populations) == 2
            [cvpairwisepopulation]
        elseif length(populations) > 2
            [cvpairwisepopulation, cvleaveonepopulationout]
        else
            throw(ErrorException("The input population in GBInput is empty."))
        end
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

Extract allele effects by fitting genomic prediction models without cross-validation.

# Arguments
- `input::GBInput`: A GBInput struct containing:
  - `fname_geno`: Path to genotype data file
  - `fname_pheno`: Path to phenotype data file
  - `models`: Vector of model functions to fit (e.g., [bayesa, bayesb])
  - `traits`: Optional vector of trait names to analyze
  - `populations`: Optional vector of population names to analyze
  - `n_iter`: Number of iterations for Bayesian models
  - `n_burnin`: Number of burn-in iterations for Bayesian models
  - `verbose`: Boolean for detailed output

# Returns
- `Vector{String}`: Paths to JLD2 files containing fitted model results, one file per model-trait-population combination

# Details
The function fits specified genomic prediction models to the full dataset without cross-validation. 
For each combination of model, trait, and population, it:
1. Loads and processes genotype and phenotype data
2. Fits the specified model
3. Saves results to a JLD2 file with naming pattern: `{prefix}_model_{name}-trait_{name}-population_{name}.jld2`

# Notes
- If populations is not specified, all unique populations in phenotype data are used
- If traits is not specified, all unique traits in phenotype data are used
- For Bayesian models (names containing "bayes"), uses specified iterations and burn-in
- Will throw an error if output files already exist

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

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
    # genomes = GenomicBreedingCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
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
    elseif !input.bulk_cv
        sort(unique(phenomes.populations))
    else
        [nothing]
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
                fname_jld2 = if isnothing(population)
                    string(fname_out_prefix, "model_", model, "-trait_", trait, "-population_bulk.jld2")
                else
                    string(fname_out_prefix, "model_", model, "-trait_", trait, "-population_", population, ".jld2")
                end
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
        
        
        
        # TODO: output this correlation matrix and/or heatmap plot




    end
    # Output filnames of the tab-delimited files containing the tables of allele effects
    fname_allele_effects_jld2s
end

"""
    predict(input::GBInput)::String

Predict trait values (GEBVs) for a set of genotypes using pre-trained models from `GenomicBreeding.fit()`.

# Arguments
- `input::GBInput`: Input configuration containing:
  - `fname_geno`: Path to genotype data file
  - `fname_allele_effects_jld2s`: Vector of paths to saved model files from previous `fit()` calls
  - `fname_out_prefix`: Prefix for output files
  - `analysis`: Set to `GenomicBreeding.predict`

# Returns
- `String`: Path to the JLD2 file containing predicted phenotypes

# Details
Takes a `GBInput` object with genotype data and pre-trained models to predict trait values 
for new individuals. The function:
1. Loads genotype data and model parameters
2. Predicts trait values using the loaded models
3. Saves predictions in a Phenomes struct
4. Returns the path to the saved predictions file

# Output File Format
- For single trait prediction: `{prefix}{model_name}-predicted_phenomes.jld2`
- For multi-trait prediction: `{prefix}{hash}-predicted_phenomes.jld2`

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

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
    # genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
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
        # println(j)
        y_pred = GenomicBreedingModels.predict(fit = model_fit, genomes = genomes, idx_entries = collect(1:n))
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

"""
    fitandpredict(input::GBInput)::Vector{String}

Performs genomic prediction by first fitting the model and then making predictions.

This function is a convenience wrapper that combines the `fit` and `predict` steps
into a single operation. It first fits the model using the provided input data,
saves the allele effects, and then generates predictions.

# Arguments
- `input::GBInput`: A GBInput object containing all necessary parameters and data for genomic prediction

# Returns
- `Vector{String}`: A vector of strings containing the prediction results

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=GenomicBreeding.fitandpredict, fname_out_prefix="GBOutput_fitandpredict/output-", verbose=false);

julia> fname_phenomes_predicted = GenomicBreeding.fitandpredict(input);

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
function fitandpredict(input::GBInput)::String
    input.fname_allele_effects_jld2s = fit(input)
    input.traits = nothing # predict all traits from the fitted models
    input.fname_out_prefix = "$(input.fname_out_prefix)-predict-"
    input.analysis = predict
    predict(input)
end