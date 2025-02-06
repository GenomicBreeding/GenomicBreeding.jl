"""
    assess(input::GBInput)::Tuple{DataFrame,DataFrame}

Assess genomic prediction accuracy via replicated k-fold cross-validation

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = writedelimited(genomes, fname="test-geno.tsv");

julia> fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv");

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], models=[ridge, bayesa], n_replications=2, n_folds=3, verbose=false);

julia> df_across, df_per_entry = assess(input);

julia> df_agg = combine(groupby(df_across, :model), :cor_mean => mean);

julia> df_agg.cor_mean_mean[df_agg.model .== "ridge"] < df_agg.cor_mean_mean[df_agg.model .== "bayesa"]
true
```
"""
function assess(input::GBInput)::Tuple{DataFrame,DataFrame}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = writedelimited(genomes, fname="test-geno.tsv")
    # fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv")
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input
    bulk_cv = input.bulk_cv
    models = input.models
    n_folds = input.n_folds
    n_replications = input.n_replications
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = load(input)
    dim_genomes = dimensions(genomes)
    dim_phenomes = dimensions(phenomes)
    # Bulk CV or if there is only one population
    cvs_bulk, notes_bulk = if bulk_cv || (dim_genomes["n_populations"] == 1)
        cvbulk(
            genomes = genomes,
            phenomes = phenomes,
            models = models,
            n_replications = n_replications,
            n_folds = n_folds,
            verbose = verbose,
        )
    else
        nothing, nothing
    end
    # Per population CV, i.e. if not builk CV (across all populations) and there are more than 1 population
    cvs_perpop, notes_perpop = if !bulk_cv && (dim_genomes["n_populations"] > 1)
        cvperpopulation(
            genomes = genomes,
            phenomes = phenomes,
            models = models,
            n_replications = n_replications,
            n_folds = n_folds,
            verbose = verbose,
        )
    else
        (nothing, nothing)
    end
    # Pairwise population CV, i.e. if not builk CV (across all populations) and there are more than 1 population
    cvs_pairwise, notes_pairwise = if !bulk_cv && (dim_genomes["n_populations"] > 1)
        cvpairwisepopulation(
            genomes = genomes,
            phenomes = phenomes,
            models = models,
            n_replications = n_replications,
            n_folds = n_folds,
            verbose = verbose,
        )
    else
        nothing, nothing
    end
    # Leave-one-population-out CV, i.e. if not builk CV (across all populations) and there are more than 1 population
    cvs_lopo, notes_lopo = if !bulk_cv && (dim_genomes["n_populations"] > 1)
        cvleaveonepopulationout(
            genomes = genomes,
            phenomes = phenomes,
            models = models,
            n_replications = n_replications,
            n_folds = n_folds,
            verbose = verbose,
        )
    else
        nothing, nothing
    end
    # Summarise
    df_across, df_per_entry = nothing, nothing
    if !isnothing(cvs_bulk)
        df_1, df_2 = summarise(cvs_bulk)
        df_1.cv_type .= "bulk"
        df_2.cv_type .= "bulk"
        if isnothing(df_across)
            df_across = df_1
            df_per_entry = df_2
        else
            df_across = [df_across; df_1]
            df_per_entry = [df_per_entry; df_2]
        end
        if verbose && length(notes_bulk) > 0
            println("Notes on bulk CV:")
            @show notes_bulk
        end
    end
    if !isnothing(cvs_perpop)
        df_1, df_2 = summarise(cvs_perpop)
        df_1.cv_type .= "perpop"
        df_2.cv_type .= "perpop"
        if isnothing(df_across)
            df_across = df_1
            df_per_entry = df_2
        else
            df_across = [df_across; df_1]
            df_per_entry = [df_per_entry; df_2]
        end
        if verbose && length(notes_perpop) > 0
            println("Notes on perpop CV:")
            @show notes_perpop
        end
    end
    if !isnothing(cvs_pairwise)
        df_1, df_2 = summarise(cvs_pairwise)
        df_1.cv_type .= "pairwise"
        df_2.cv_type .= "pairwise"
        if isnothing(df_across)
            df_across = df_1
            df_per_entry = df_2
        else
            df_across = [df_across; df_1]
            df_per_entry = [df_per_entry; df_2]
        end
        if verbose && length(notes_pairwise) > 0
            println("Notes on pairwise CV:")
            @show notes_pairwise
        end
    end
    if !isnothing(cvs_lopo)
        df_1, df_2 = summarise(cvs_lopo)
        df_1.cv_type .= "lopo"
        df_2.cv_type .= "lopo"
        if isnothing(df_across)
            df_across = df_1
            df_per_entry = df_2
        else
            df_across = [df_across; df_1]
            df_per_entry = [df_per_entry; df_2]
        end
        if verbose && length(notes_lopo) > 0
            println("Notes on lopo CV:")
            @show notes_lopo
        end
    end
    # Output
    (df_across, df_per_entry)
end

"""
    extracteffects(input::GBInput)::Vector{DataFrame}

Extract allele effects by fitting the models without cross-validation

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = writedelimited(genomes, fname="test-geno.tsv");

julia> fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv");

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], models=[ridge, bayesa], n_replications=2, n_folds=3, verbose=false);

julia> effects = extracteffects(input);

julia> sum([nrow(effects[i]) == nrow(effects[1]) for i in eachindex(effects)]) == length(effects)
true
```
"""
function extracteffects(input::GBInput)::Vector{DataFrame}
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
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = load(input)
    # Fit the entire data to extract effects per trait per model
    populations = sort(unique(phenomes.populations))
    traits = phenomes.traits
    out::Vector{DataFrames.DataFrame} = fill(DataFrame(), (1 + length(populations)) * length(traits) * length(models))
    if verbose
        pb = ProgressMeter.Progress(length(out); desc="Extracting allele effects: ")
    end
    for (i, model) in enumerate(models)
        for (j, trait) in enumerate(traits)
            for (k, population) in enumerate(vcat("bulk", populations))
                # i = 1; model = models[i]; j = 1; trait = traits[j]; k = 1; population = populations[k]
                idx = ((((i-1)*length(traits) + j) - 1) * (1+length(populations))) + k
                println(idx)
                if verbose
                    println(string("Model: ", model, "| Trait: ", trait, "| Population: ", population))
                end
                Γ, Φ = if population == "bulk"
                    (genomes, phenomes)
                else
                    (
                        slice(genomes, idx_entries=findall(genomes.populations .== population)),
                        slice(phenomes, idx_entries=findall(phenomes.populations .== population)),
                    )
                end
                fit = if !isnothing(match(Regex("bayes", "i"), string(model)))
                    model(genomes=Γ, phenomes=Φ, n_iter=n_iter, n_burnin=n_burnin, verbose=false)
                else
                    model(genomes=Γ, phenomes=Φ, verbose=false)
                end
                out[idx] = DataFrame(
                    model=fit.model,
                    trait=fit.trait,
                    population=population,
                    fit_corr=fit.metrics["cor"],
                    b_hat_labels=fit.b_hat_labels,
                    b_hat=fit.b_hat,
                )
                if verbose
                    ProgressMeter.next!(pb)
                end
            end
        end
    end
    if verbose
        ProgressMeter.finish!(pb)
        ### Correlations between in allele effects
        C = fill(-Inf, length(out), length(out))
        for i in eachindex(out)
            for j in eachindex(out)
                C[i, j] = cor(out[i].b_hat[2:end], out[j].b_hat[2:end])
                # model1 = out[i].model[1]; trait1 = out[i].trait[1]; population1 = out[i].population[1]
                # model2 = out[j].model[1]; trait2 = out[j].trait[1]; population2 = out[j].population[1]
                # println(join([model1, trait1, population1], "|"), " vs ",join([model2, trait2, population2], "|"), ": ", round(C[i, j], digits=2))
            end
        end
        p = UnicodePlots.heatmap(C, title="Correlation between allele effects")
        display(p)
    end
    # Output
    out
end


