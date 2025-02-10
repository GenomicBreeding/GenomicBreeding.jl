function plot(input::GBInput)::String
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # # fname_geno = "test-geno.tsv"; fname_pheno = "test-pheno.tsv"; outdir = "GBOutput"
    # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # # input.fname_out_prefix = replace(input.fname_out_prefix, dirname(input.fname_out_prefix) => outdir)
    # Load genomes and phenomes
    genomes, phenomes = load(input)
    # Load CVs
    directory_name = if dirname(input.fname_out_prefix) == ""
        pwd()
    else
        dirname(input.fname_out_prefix)
    end
    files = readdir(directory_name)
    idx = findall(.!isnothing.(match.(Regex("jld2\$"), files)))
    fnames_cvs = if length(idx) > 0
        joinpath.(directory_name, files[idx])
    else
        nothing
    end
    
    # Plot phenotypes
    distribution_plots = GBPlots.plot(DistributionPlots, phenomes)
    violin_plots = GBPlots.plot(ViolinPlots, phenomes)
    corheat_plots = GBPlots.plot(CorHeatPlots, phenomes)
    tree_plots = GBPlots.plot(TreePlots, phenomes)
    
    
    # Plot genotypes
    
    
    # Plot CVs, if present
    bar_plots = if !isnothing(fnames_cvs)
        cvs = Vector{CV}(undef, length(fnames_cvs))
        for (i, fname) in enumerate(fnames_cvs)
            # i = 1; fname = fnames_cvs[i];
            cvs[i] = readjld2(CV, fname=fname)
        end
        GBPlots.plot(BarPlots, cvs)
    else
        nothing
    end
    bar_plots.plots[1]

    # Save plots



end
