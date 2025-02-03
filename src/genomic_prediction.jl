function assess(input::GBInput)::Tuple{DataFrame,DataFrame}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = writedelimited(genomes, fname="test-geno.tsv")
    # fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv")
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input
    # fname_geno = input.fname_geno
    # fname_pheno = input.fname_pheno
    bulk_cv = input.bulk_cv
    # populations = input.populations
    # traits = input.traits
    models = input.models
    n_folds = input.n_folds
    n_replications = input.n_replications
    # keep_all = input.keep_all
    # maf = input.maf
    # mtv = input.mtv
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
