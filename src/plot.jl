function plot(input::GBInput)::String
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # # directory_name = submitslurmarrayjobs(input=input, analysis=assess)
    # Load genomes and phenomes
    genomes, phenomes = load(input)


    # Plot phenotypes
    GBPlots.plot(DistributionPlots, phenomes)
    GBPlots.plot(ViolinPlots, phenomes)
    GBPlots.plot(CorHeatPlots, phenomes)
    GBPlots.plot(TreePlots, phenomes)
    # Plot genotypes
    # Plot CVs, if present
end