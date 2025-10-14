"""
    submitslurmarrayjobs(; input::GBInput, auto_proceed::Bool=false)::String

Submit an array of Slurm jobs for genomic prediction analysis.

# Arguments
- `input::GBInput`: A GBInput struct containing all necessary parameters for job submission and analysis.
- `auto_proceed::Bool=false`: If true, skips the interactive confirmation prompt before job submission.

# Returns
- `String`: Path to the output directory where results will be stored.

# Details
This function handles the submission of parallel genomic prediction jobs to a Slurm cluster. It performs the following steps:
1. Validates input parameters and checks for required R packages
2. Creates necessary output directories
3. Prepares individual job inputs
4. Generates Julia and Slurm scripts
5. Submits the array job to the Slurm scheduler

The function supports various genomic analyses including:
- Cross-validation (cv)
- Model fitting (fit)
- Prediction (predict)
- GWAS analysis (gwas)

# Job Configuration
- Uses Slurm array jobs for parallel execution
- Configurable CPU, memory, and time limit parameters
- Supports both module-based and conda environments
- Interactive confirmation before job submission

# Notes
- Requires a working Slurm environment
- BGLR R package must be installed
- User will be prompted to enter "YES" to confirm job submission
- Job array size is controlled by `SLURM_max_array_jobs_running`

# Example
```julia
using GenomicBreeding, StatsBase;
using GenomicBreeding: cv, fit, predict, fitandpredict, gwas, ols, rigde, lasso, bayesa, bayesb, bayesc, gwasols, gwaslmm, gwasreml;
genomes = GenomicBreeding.GenomicBreedingCore.simulategenomes(n=300, l=1_000, n_populations=3, verbose=false);
trials, _ = GenomicBreeding.GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
phenomes = extractphenomes(trials);
fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

# Repeated k-fold cross-validation
input_cv = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, SLURM_account_name="dbiof1", SLURM_cpus_per_task=5, SLURM_mem_G=5, fname_out_prefix="GBOutput/test-", verbose=false);
outdir = submitslurmarrayjobs(input_cv); ### You will be asked to enter "YES" to proceed with job submission.
run(`sh -c 'squeue -u "\$USER"'`)
run(`sh -c 'tail slurm-*_*.out'`)
run(`sh -c 'grep -i "err" slurm-*_*.out'`)
cvs = loadcvs(input_cv)
df_across_entries, df_per_entry = tabularise(cvs)
sum(df_across_entries.model .== "bayesa") / nrow(df_across_entries)
sum(df_across_entries.model .== "ridge") / nrow(df_across_entries)
sort(combine(groupby(df_across_entries, [:validation_population, :model]), [:cor => mean, :cor => length]), :cor_mean, rev=true)

# Genomic prediction equation full data fit
input_fit = clone(input_cv)
input_fit.analysis = fit
outdir = submitslurmarrayjobs(input_fit);
input_fit.fname_allele_effects_jld2s = begin
    files = readdir(outdir)
    idx = findall(.!isnothing.(match.(Regex("-fit-"), files)) .&& .!isnothing.(match.(Regex("jld2\$"), files)))
    joinpath.(outdir, files[idx])
end
fits = loadfits(input_fit)
length(fits)

# Calculate GEBVs
input_predict = clone(input_fit)
input_predict.analysis = predict
outdir = submitslurmarrayjobs(input_predict)
run(`squeue`)

# Fit and predict in one step
input_fitandpredict = clone(input_fit)
input_fitandpredict.analysis = fitandpredict
outdir = submitslurmarrayjobs(input_fitandpredict, auto_proceed=true)
run(`squeue`)
```
"""
function submitslurmarrayjobs(input::GBInput; auto_proceed::Bool=false)::String
    # genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv, SLURM_cpus_per_task=6, SLURM_mem_G=9)
    # Catch the main errors early
    errors::Vector{String} = checkinputs(input)
    begin
        commands = if input.SLURM_module_load_R_version_name != "conda"
            [
                "#!/bin/bash",
                string("module load ", input.SLURM_module_load_R_version_name),
                "Rscript -e 'if(sum(rownames(installed.packages()) == \"BGLR\") > 0){cat(0)}else{cat(404)}'",
            ]
        else
            [
                "#!/bin/bash",
                string("module load ", input.SLURM_module_load_Conda_version_name),
                "conda init bash",
                "source ~/.bashrc",
                "conda activate GenomicBreeding",
                "Rscript -e 'if(sum(rownames(installed.packages()) == \"BGLR\") > 0){cat(0)}else{cat(404)}'",
            ]
        end
        open("test-check_BGLR.sh", "w") do file
            write(file, join(commands, "\n") * "\n")
        end
        run(`chmod +x test-check_BGLR.sh`)
        out = Pipe()
        try
            run(pipeline(`./test-check_BGLR.sh`, stdout = out))
        catch
            throw(
                ErrorException(
                    "Please specify the correct R module. Invoked the followinfg with an error: `module load " *
                    input.SLURM_module_load_R_version_name *
                    "`.",
                ),
            )
        end
        close(out.in)
        rm("test-check_BGLR.sh")
        if String(read(out)) == "404"
            push!(
                errors,
                string(
                    "The BGLR R package is not installed. Please install it manually first (e.g. module load R; Rscript -e 'install.packages(\"BGLR\")')",
                ),
            )
        end
    end
    if length(errors) > 0
        throw(ArgumentError("Error/s:\n\t‣ " * join(errors, "\n\t‣ ")))
    end
    # Define the prefix of the output including the output directory
    fname_out_prefix = input.fname_out_prefix
    # Instantiate the output directory if it does not exist
    outdir = if dirname(fname_out_prefix) == ""
        pwd()
    else
        dirname(fname_out_prefix)
    end
    if !isdir(outdir)
        try
            mkdir(outdir)
        catch
            throw(ArgumentError("Error creating the output directory: `" * outdir * "`"))
        end
    end
    # Create a run directory in the directory defined in the fname_out_prefix
    run_outdir = joinpath(outdir, "run")
    if !isdir(run_outdir)
        try
            mkdir(run_outdir)
        catch
            throw(ArgumentError("Error creating the output directory: `" * run_outdir * "`"))
        end
    end
    # Check the input files, define the vector of GBInput structs for parallel execution and save them in the run directory
    inputs = prepareinputs(input)
    n_array_jobs = length(inputs)
    for i = 1:n_array_jobs
        # i = 1
        writejld2(inputs[i], fname = joinpath(run_outdir, string("GBInput-", input.analysis, "-", i, ".jld2")))
    end
    # Save the Julia run file in the run directory
    julia_script = [
        "i = ARGS[1]",
        # "try",
        # "\tusing Pkg; Pkg.update()",
        # "catch",
        # "\tnothing",
        # "end",
        "using GenomicBreeding",
        "import GenomicBreeding: cv, fit, predict, fitandpredict, gwas",
        "import GenomicBreeding: ols, ridge, lasso, bayesa, bayesb, bayesc, gwasols, gwaslmm, gwasreml",
        string(
            "input = readjld2(GBInput, fname=joinpath(\"",
            run_outdir,
            "\", string(\"GBInput-\", ",
            input.analysis,
            ", \"-\", i, \".jld2\")))",
        ),
        "display(input)",
        string("output = ", input.analysis, "(input)"),
        "display(output)",
    ]
    open(joinpath(run_outdir, "run.jl"), "w") do file
        write(file, join(julia_script, "\n") * "\n")
    end
    # Define the Slurm run file
    slurm_script = [
        "#!/bin/bash",
        string("#SBATCH --job-name='", input.SLURM_job_name, "'"),
        string("#SBATCH --account='", input.SLURM_account_name, "'"),
        string("#SBATCH --partition='", input.SLURM_partition_name, "'"),
        string("#SBATCH --nodes=", input.SLURM_nodes_per_array_job),
        string("#SBATCH --ntasks=", input.SLURM_tasks_per_node),
        string("#SBATCH --cpus-per-task=", input.SLURM_cpus_per_task),
        string("#SBATCH --gres=gpu:", input.SLURM_gpus),
        string("#SBATCH --mem=", input.SLURM_mem_G, "G"),
        string("#SBATCH --time=", input.SLURM_time_limit_dd_hhmmss),
        "LD_LIBRARY_PATH=\"\"",
        "export JULIA_CPU_TARGET=generic", # to prevent constant precompilations
        string("module load ", input.SLURM_module_load_R_version_name),
        string("module load ", input.SLURM_module_load_Julia_version_name),
        string(
            "time julia --threads ",
            input.SLURM_cpus_per_task,
            " ",
            joinpath(run_outdir, "run.jl \$SLURM_ARRAY_TASK_ID"),
        ),
    ]
    # Remove the account name, partion and Julia module (use manually installed Julia) if empty
    if input.SLURM_account_name == ""
        slurm_script = slurm_script[isnothing.(match.(Regex("^#SBATCH --account="), slurm_script))]
    end
    if input.SLURM_gpus == 0
        slurm_script = slurm_script[isnothing.(match.(Regex("^#SBATCH --gres=gpu:"), slurm_script))]
    end
    if input.SLURM_partition_name == ""
        slurm_script = slurm_script[isnothing.(match.(Regex("^#SBATCH --partition="), slurm_script))]
    end
    if input.SLURM_module_load_R_version_name == "conda"
        idx = findall(.!isnothing.(match.(Regex("^module load conda\$"), slurm_script)))[1]
        slurm_script[idx] = string(
            "module load ",
            input.SLURM_module_load_Conda_version_name,
            ";",
            "conda init bash;",
            "source ~/.bashrc;",
            "conda activate GenomicBreeding;",
        )
    end
    if input.SLURM_module_load_Julia_version_name == ""
        slurm_script = slurm_script[isnothing.(match.(Regex("^module load \$"), slurm_script))]
    end
    # Save the Slurm run file
    open(joinpath(run_outdir, "run.slurm"), "w") do file
        write(file, join(slurm_script, "\n") * "\n")
    end
    if !auto_proceed
        # Ask the use interactively to confirm
        println(
            string(
                "Are you sure you want to submit a total of ",
                n_array_jobs,
                " jobs (with a maximum of ",
                input.SLURM_max_array_jobs_running,
                " jobs running simultaneously) each requiring ",
                input.SLURM_cpus_per_task,
                " cpus and ",
                input.SLURM_mem_G,
                "Gb of RAM with a time limit of ",
                input.SLURM_time_limit_dd_hhmmss,
                "?",
            ),
        )
        println("Please enter YES to proceed, otherwise enter anything else to cancel:")
        proceed = strip(readline())
        @show proceed == "YES"
        if proceed != "YES"
            println(proceed)
            println("Cancelled.")
            println("If this was a mistake you may submit the array jobs manually via:")
            println(
                join(
                    [
                        "sbatch",
                        string("--array=1-", n_array_jobs, "%", input.SLURM_max_array_jobs_running),
                        joinpath(run_outdir, "run.slurm"),
                    ],
                    " ",
                ),
            )
            return outdir
        end
    end
    # Submit the array of jobs
    SHELL_COMMAND = Cmd([
        "sbatch",
        string("--array=1-", n_array_jobs, "%", input.SLURM_max_array_jobs_running),
        joinpath(run_outdir, "run.slurm"),
    ])
    run(SHELL_COMMAND)
    # Output the name of the output directory
    outdir
end
