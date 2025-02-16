var documenterSearchIndex = {"docs":
[{"location":"models/","page":"Models","title":"Models","text":"CurrentModule = GenomicBreeding","category":"page"},{"location":"models/#Models","page":"Models","title":"Models","text":"","category":"section"},{"location":"models/","page":"Models","title":"Models","text":"Models","category":"page"},{"location":"models/","page":"Models","title":"Models","text":"","category":"page"},{"location":"models/#Linear-mixed-models-for-analysing-trial","page":"Models","title":"Linear mixed models for analysing trial","text":"","category":"section"},{"location":"models/#Linear-models-for-genomic-prediction","page":"Models","title":"Linear models for genomic prediction","text":"","category":"section"},{"location":"models/#Non-linear-models-for-genomic-prediction","page":"Models","title":"Non-linear models for genomic prediction","text":"","category":"section"},{"location":"models/#Epifeat","page":"Models","title":"Epifeat","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GenomicBreeding","category":"page"},{"location":"#GenomicBreeding","page":"Home","title":"GenomicBreeding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GenomicBreeding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GenomicBreeding]","category":"page"},{"location":"#GenomicBreeding.GBInput","page":"Home","title":"GenomicBreeding.GBInput","text":"mutable struct GBInput\n    fname_geno::String\n    fname_pheno::String\n    fname_allele_effects_jld2s::Vector{String}\n    analysis::Function\n    bulk_cv::Bool\n    populations::Union{Nothing,Vector{String}}\n    traits::Union{Nothing,Vector{String}}\n    models::Any\n    n_folds::Int64\n    n_replications::Int64\n    gwas_models::Any\n    keep_all::Bool\n    maf::Float64\n    mtv::Float64\n    n_iter::Int64\n    n_burnin::Int64\n    fname_out_prefix::String\n    SLURM_job_name::String\n    SLURM_account_name::String\n    SLURM_partition_name::String\n    SLURM_nodes_per_array_job::Int64\n    SLURM_tasks_per_node::Int64\n    SLURM_cpus_per_task::Int64\n    SLURM_mem_G::Float64\n    SLURM_time_limit_dd_hhmmss::String\n    SLURM_max_array_jobs_running::Int64\n    SLURM_module_load_R_version_name::String\n    SLURM_module_load_BLAS_version_name::String\n    verbose::Bool\nend\n\nInput struct (belongs to GBCore.AbstractGB type)\n\nfname_geno: genotype file (see file format guide: TODO: {URL HERE})\nfname_pheno: phenotype file (see file format guide: TODO: {URL HERE}; Default = \"\")\nfname_allele_effects_jld2s: vector of filenames of JLD2 files containing the Fit struct of a genomic prediction model (Default = [\"\"])\nanalysis: analysis to perform or function to use (Default = cv):\ncv: replicated k-fold cross-validation\nfit: fit genomic prediction models without cross-validation to extract allele effects to compute GEBVs on other genomes\npredict: compute GEBVs using the output of fit and genotype data lacking empirical GP-model-associated phenotype data in the fit output (do not forget to set the fname_allele_effects_jld2s field)\ngwas: genome-wide association study\nbulk_cv: perform cross-validation across all populations, i.e. disregard population grouping (Default = false)\npopulations: include only these populations (Default = nothing which means include all populations)\ntraits: include only these traits (Default = nothing which means include all traits)\nmodels: genomic prediction model functions (Default = [ridge, bayesa]; see models list: TODO: {URL HERE})\nn_folds: number of k partitions for k-fold cross-validation (Default = 5)\nn_replications: number of replications for repeated k-fold cross-validation (Default = 5)\ngwas_models: models to use if performning genome-wide association study (Default = [gwasols, gwaslmm])\nkeep_all: keep all entries upon merging genomes and phenomes potentially resulting in sparsities in both structs? (Default = false)\nmaf: minimum allele frequency (Default = 0.05)\nmtv: minimum trait variance (Default = 1e-7)\nn_iter: number of Bayesian model fitting MCMC/HMC iteration (Default = 1_500)\nn_burnin: number of initial Bayesian model fitting MCMC/HMC iterations to be excluded from the posterior distribution (Default = 500)\nfname_out_prefix: prefix of the output files which may include directory names (Default = \"\" which translates to ./GBOuput/output-<yyyymmddHHMMSS>-<3_digit_random_number>-)\nSLURM_job_name: name of the Slurm job array (Default = \"\")\nSLURM_account_name: Slurm account name (Default = \"\")\nSLURM_partition_name: Slurm paritition to use (Default = \"\")\nSLURM_nodes_per_array_job: number of nodes per array job (Default = 1)\nSLURM_tasks_per_node: number of tasks per node (Default = 1)\nSLURM_cpus_per_task: number of CPU cores per task (Default = 1)\nSLURM_mem_G: maximum memroy requested in gigabytes (Default = 1.00)\nSLURM_time_limit_dd_hhmmss: maximum computation time requested in the follwowing format: \"dd-hh:mm:ss\" (Default = \"00-01:00:00\")\nSLURM_max_array_jobs_running: maximum number of array jobs which can run simultaneously (Default = 20)\nSLURM_module_load_R_version_name: name of the R statistical language module and version to be used to call R::BGLR (Default = \"R\")\nSLURM_module_load_BLAS_version_name: name of the BLAS module and version to be used to use with R::BGLR. This is non-critical as there is usually a system default to fallback to. (Default = \"\")\nverbose: show messages (Default = true)\n\n\n\n\n\n","category":"type"},{"location":"#Base.:==-Tuple{GBInput, GBInput}","page":"Home","title":"Base.:==","text":"Base.:(==)(x::GBInput, y::GBInput)::Bool\n\nEquality of GBInput structs using the hash function defined for GBInput structs.\n\nExamples\n\njulia> input_1 = input = GBInput(fname_geno=\"geno1.jld2\", fname_pheno=\"pheno1.jld2\", fname_out_prefix=\"test1-\", SLURM_job_name=\"slurmjob1\");\n\njulia> input_2 = input = GBInput(fname_geno=\"geno1.jld2\", fname_pheno=\"pheno1.jld2\", fname_out_prefix=\"test1-\", SLURM_job_name=\"slurmjob1\");\n\njulia> input_3 = input = GBInput(fname_geno=\"geno2.jld2\", fname_pheno=\"pheno2.jld2\");\n\njulia> input_1 == input_2\ntrue\n\njulia> input_1 == input_3\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#Base.hash-Tuple{GBInput, UInt64}","page":"Home","title":"Base.hash","text":"Base.hash(x::GBInput, h::UInt)::UInt\n\nHash a GBInput struct using all its fields. We deliberately excluded the allele_frequencies, and mask for efficiency.\n\nExamples\n\njulia> input = GBInput(fname_geno=\"\", fname_pheno=\"\");\n\njulia> typeof(hash(input))\nUInt64\n\n\n\n\n\n","category":"method"},{"location":"#GBCore.checkdims-Tuple{GBInput}","page":"Home","title":"GBCore.checkdims","text":"checkdims(input::GBInput)::Bool\n\nCheck dimension compatibility of the fields of the GBInput struct\n\nExamples\n\njulia> input = GBInput(fname_geno=\"geno1.jld2\", fname_pheno=\"pheno1.jld2\");\n\njulia> checkdims(input)\ntrue\n\njulia> input.models = nothing\n\njulia> checkdims(input)\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.checkinputs-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.checkinputs","text":"checkinputs(input::GBInput)::Bool\n\nCheck input compatibility with the analysis requested\n\nExamples\n\njulia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n    \njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv);\n\njulia> length(checkinputs(input)) == 0\ntrue\n\njulia> input.fname_pheno = \"\"; length(checkinputs(input)) == 0\nfalse\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.clone-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.clone","text":"clone(x::GBInput)::GBInput\n\nClone a GBInput object\n\nExample\n\njulia> input = GBInput(fname_geno=\"geno1.jld2\", fname_pheno=\"pheno1.jld2\");\n\njulia> copy_input = clone(input);\n\njulia> input == copy_input\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.cv-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.cv","text":"cv(input::GBInput)::Tuple{Vector{String},Vector{String}}\n\nAssess genomic prediction accuracy via replicated k-fold cross-validation. Outputs are saved as JLD2 (each containing a CV struct per fold, replication, and trait) and possibly text file/s containing notes describing why some jobs failed.\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=[\"pop_1\", \"pop_3\"], traits=[\"trait_1\"], n_replications=2, n_folds=3, verbose=false);\n\njulia> fnames_cvs, fnames_notes = cv(input);\n\njulia> length(fnames_cvs) == 32, length(fnames_notes) == 0\n(true, true)\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.fit-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.fit","text":"fit(input::GBInput)::Vector{String}\n\nExtract allele effects by fitting the models without cross-validation. Outputs are JLD2 files, each containing the Fit struct for each model-trait-population combination.\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=[\"pop_1\", \"pop_3\"], traits=[\"trait_1\"], n_replications=2, n_folds=3, verbose=false);\n\njulia> fname_allele_effects_jld2s = GenomicBreeding.fit(input);\n\njulia> length(fname_allele_effects_jld2s) == 6\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.gwas-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.gwas","text":"gwas(input::GBInput)::Vector{String}\n\nPerform genome-wide association study. Outputs are JLD2 files, each containing the Fit struct for each model-trait-population combination. Note that the b_hat field of each Fit struct contains t-statistics for gwasols and z-statistics for the other GWAS models, instead of allele effects in genomic prediction models.\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, gwas_models=[gwasols], verbose=false);\n\njulia> fname_test_statistics_jld2s = GenomicBreeding.gwas(input);\n\njulia> length(fname_test_statistics_jld2s) == 12\ntrue\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.load-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.load","text":"load(input::GBInput)::Tuple{Genomes, Phenomes}\n\nLoad, merge and filter genotype and phenotype data\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n    \njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=[\"pop_1\", \"pop_3\"], traits=[\"trait_1\"], verbose=false);\n\njulia> genomes, phenomes = load(input);\n\njulia> length(unique(genomes.populations)) == length(unique(phenomes.populations)) == 2\ntrue\n\njulia> length(phenomes.traits) == 1\ntrue\n\njulia> rm.([fname_geno, fname_pheno]);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.plot-Tuple{}","page":"Home","title":"GenomicBreeding.plot","text":"plot(;\n    input::GBInput,\n    skip_genomes::Bool = false,\n    skip_phenomes::Bool = false,\n    format::String = \"svg\",\n    plot_size::Tuple{Int64,Int64} = (600, 450),\n    overwrite::Bool = false,\n)::String\n\nPlot genomes, phenomes, and CVs, if present.\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n    \njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5, fname_out_prefix=\"GBOutput/test-\", verbose=false);\n\njulia> GenomicBreeding.plot(input=input, format=\"png\", plot_size = (700, 525))\n\"GBOutput/plots\"\n\njulia> GenomicBreeding.plot(input=input, format=\"png\", plot_size = (700, 525), overwrite=true, skip_genomes=true)\n\"GBOutput/plots\"\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.predict-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.predict","text":"fit(input::GBInput)::Vector{String}\n\nExtract allele effects by fitting the models without cross-validation. Generated JLD2 files each containing the Fit struct for each model-trait-population combination.\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, verbose=false);\n\njulia> input.fname_allele_effects_jld2s = GenomicBreeding.fit(input);\n\njulia> fname_phenomes_predicted = GenomicBreeding.predict(input);\n\njulia> phenomes_predicted = readjld2(Phenomes, fname=fname_phenomes_predicted);\n\njulia> dimensions(phenomes_predicted)\nDict{String, Int64} with 8 entries:\n  \"n_total\"       => 7200\n  \"n_zeroes\"      => 0\n  \"n_nan\"         => 0\n  \"n_entries\"     => 300\n  \"n_traits\"      => 24\n  \"n_inf\"         => 0\n  \"n_populations\" => 3\n  \"n_missing\"     => 0\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.prepareinputs-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.prepareinputs","text":"prepareinputs(input::GBInput)::Vector{GBInput}\n\nPrepare GBInputs for Slurm array jobs\n\nExample\n\njulia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n    \njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, verbose=false);\n\njulia> inputs = prepareinputs(input);\n\njulia> length(inputs) == 30\ntrue\n\njulia> rm.([fname_geno, fname_pheno]);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.prepareoutprefixandoutdir-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.prepareoutprefixandoutdir","text":"prepareoutprefixandoutdir(input::GBInput)::String\n\nPrepare the output prefix by replacing problematic strings in the prefix of the output filenames and instantiating the output folder\n\nExample\n\njulia> input = GBInput(fname_geno=\"some_dir/fname_geno.jld2\", fname_pheno=\"some_dir/fname_pheno.jld2\", fname_out_prefix=\"GBOutput/some@!_%&prefix-\", verbose=false);\n\njulia> fname_out_prefix = prepareoutprefixandoutdir(input)\n\"GBOutput/some_____prefix-\"\n\njulia> rm(dirname(fname_out_prefix), recursive=true);\n\n\n\n\n\n","category":"method"},{"location":"#GenomicBreeding.submitslurmarrayjobs-Tuple{GBInput}","page":"Home","title":"GenomicBreeding.submitslurmarrayjobs","text":"submitslurmarrayjobs(; input::GBInput)::String\n\nAssess genomic prediction accuracy via replicated k-fold cross-validation. Outputs are saved as JLD2 (each containing a CV struct per fold, replication, and trait) and possibly text file/s containing notes describing why some jobs failed. Note that you will be prompted to enter YES to proceed with Slurm job submission after showing you the job details to review and confirm.\n\nExample\n\n<!– ```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames) –>\n\njulia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.(\"pop_\", 1:3), length(genomes.entries), replace=true);\n\njulia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);\n\njulia> phenomes = extractphenomes(trials);\n\njulia> fname_geno = try writedelimited(genomes, fname=\"test-geno.tsv\"); catch; rm(\"test-geno.tsv\"); writedelimited(genomes, fname=\"test-geno.tsv\"); end;\n\njulia> fname_pheno = try writedelimited(phenomes, fname=\"test-pheno.tsv\"); catch; rm(\"test-pheno.tsv\"); writedelimited(phenomes, fname=\"test-pheno.tsv\"); end;\n\njulia> input_cv = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, SLURM_cpus_per_task=1, SLURM_mem_G=0.5, fname_out_prefix=\"GBOutput/test-\", verbose=false);\n\njulia> outdir = submitslurmarrayjobs(input_cv);\n\njulia> input_fit = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=fit, SLURM_cpus_per_task=1, SLURM_mem_G=0.5, fname_out_prefix=\"GBOutput/test-\", verbose=false);\n\njulia> outdir = submitslurmarrayjobs(input_fit);\n\n\n\n\njulia> input_predict = GBInput(fname_geno=fname_geno, fname_allele_effects_jld2s=fname_allele_effects_jld2s analysis=predict, SLURM_cpus_per_task=1, SLURM_mem_G=0.5, fname_out_prefix=\"GBOutput/test-\", verbose=false);\n\njulia> input_cv = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, SLURM_cpus_per_task=1, SLURM_mem_G=0.5, fname_out_prefix=\"GBOutput/test-\", verbose=false);\n\njulia> outdir = submitslurmarrayjobs(input_cv)\nGBOutput\n\njulia> run(`squeue`)\n\n\n\n\n\n","category":"method"},{"location":"simulations/","page":"Simulation","title":"Simulation","text":"CurrentModule = GenomicBreeding","category":"page"},{"location":"simulations/#Simulation","page":"Simulation","title":"Simulation","text":"","category":"section"},{"location":"simulations/","page":"Simulation","title":"Simulation","text":"Simulations","category":"page"},{"location":"simulations/","page":"Simulation","title":"Simulation","text":"","category":"page"},{"location":"simulations/#Simulate-genomes","page":"Simulation","title":"Simulate genomes","text":"","category":"section"},{"location":"simulations/#Simulate-effects","page":"Simulation","title":"Simulate effects","text":"","category":"section"},{"location":"simulations/#Simulate-trials","page":"Simulation","title":"Simulate trials","text":"","category":"section"}]
}
