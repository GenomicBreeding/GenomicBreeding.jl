function gwas(input::GBInput)::Vector{String}
    # genomes = GBCore.simulategenomes(n=300, l=100, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, traits=["trait_1"])
    # Parse input
    models = input.gwas_models
    fname_out_prefix = input.fname_out_prefix
    verbose = input.verbose
    # Load genomes and phenomes
    genomes, phenomes = load(input)
    dim_genomes = dimensions(genomes)
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
    # GWAS

    # function gwasols(
    #     genomes::Genomes,
    #     phenomes::Phenomes;
    #     idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_loci_alleles::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_trait::Int64 = 1,
    #     GRM_type::String = ["simple", "ploidy-aware"][1],
    #     verbose::Bool = false,
    # )::Fit
    # function gwaslmm(
    #     genomes::Genomes,
    #     phenomes::Phenomes;
    #     idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_loci_alleles::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_trait::Int64 = 1,
    #     GRM_type::String = ["simple", "ploidy-aware"][1],
    #     verbose::Bool = false,
    # )::Fit
    # function gwasreml(
    #     genomes::Genomes,
    #     phenomes::Phenomes;
    #     idx_entries::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_loci_alleles::Union{Nothing,Vector{Int64}} = nothing,
    #     idx_trait::Int64 = 1,
    #     GRM_type::String = ["simple", "ploidy-aware"][1],
    #     verbose::Bool = false,
    # )::Fit

end
