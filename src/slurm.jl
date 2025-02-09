function submitslurmarrayjobs(; input::GBInput, analysis::Function)::Nothing
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # analysis = assess
    # Catch the 3 main errors early:
    #   1. input genotype file exists
    #   2. input phenotype file exists
    #   3. R::BGLR package is installed.
    errors::Vector{String} = []
    if !isfile(input.fname_geno)
        push!(errors, string("The genotype file (", input.fname_geno, ") does not exist."))
    end
    if !isfile(input.fname_pheno)
        push!(errors, string("The phenotype file (", input.fname_pheno, ") does not exist."))
    end
    push!(
        errors,
        begin
            commands = [
                "#!/bin/bash",
                string("module load ", input.SLURM_module_load_R_version_name),
                "Rscript -e 'if(sum(rownames(installed.packages()) == \"BGLR\") > 0){cat(0)}else{cat(404)}'",
            ]
            open("test-check_BGLR.sh", "w") do file
                write(file, join(commands, "\n") * "\n")
            end
            run(`chmod +x test-check_BGLR.sh`)
            out = Pipe()
            run(pipeline(`./test-check_BGLR.sh`, stdout = out))
            close(out.in)
            rm("test-check_BGLR.sh")
            err = if String(read(out)) == "404"
                string("The BGLR R package is not installed.")
            else
                ""
            end
            err
        end,
    )
    if length(errors) > 0
        throw(ArgumentError("Error/s:\n\t‣ " * join(errors, "\n\t‣ ")))
    end
    # Define the prefix of the output including the output directory
    fname_out_prefix = input.fname_out_prefix
    # Instantiate the output directory if it does not exist
    directory_name = dirname(fname_out_prefix)
    if !isdir(directory_name) && (directory_name != "")
        try
            mkdir(directory_name)
        catch
            throw(ArgumentError("Error creating the output directory: `" * directory_name * "`"))
        end
    end
    # Create a run directory in the directory defined in the fname_out_prefix
    run_directory_name = joinpath(directory_name, "run")
    if !isdir(run_directory_name)
        try
            mkdir(run_directory_name)
        catch
            throw(ArgumentError("Error creating the output directory: `" * run_directory_name * "`"))
        end
    end
    # Check the input files, define the vector of GBInput structs for parallel execution and save them in the run directory
    inputs = prepareinputs(input)
    n_array_jobs = length(inputs)
    for i = 1:n_array_jobs
        # i = 1
        writejld2(inputs[i], fname = joinpath(run_directory_name, string("GBInput-", i, ".jld2")))
    end
    # Save the Julia run file in the run directory
    julia_script = [
        "i = ARGS[1]",
        "using GenomicBreeding",
        "import GenomicBreeding: ols, ridge, lasso, bayesa, bayesb, bayesc",
        string(
            "input = readjld2(GBInput, fname=joinpath(\"",
            run_directory_name,
            "\", string(\"GBInput-\", i, \".jld2\")))",
        ),
        "output = analysis(input)",
        "display(output)",
    ]
    open(joinpath(run_directory_name, "run.jl"), "w") do file
        write(file, join(julia_script, "\n") * "\n")
    end
    # Save the Slurm run file
    input.SLURM_job_name = string("GBJob-", Dates.format(now(), "yyyymmddHHMMSS"), Int(round(1_000 * rand())))
    input.SLURM_account_name = "some_account"
    input.SLURM_partition_name = "some_partition"
    input.SLURM_nodes_per_array_job = 1
    input.SLURM_tasks_per_node = 1
    input.SLURM_cpus_per_task = 1
    input.SLURM_mem_G = 1
    input.SLURM_time_limit_dd_hhmmss = "00-01:00:00"
    input.SLURM_max_array_jobs_running = 20
    input.SLURM_module_load_R_version_name = "R"
    slurm_script = [
        string("#SBATCH --job-name=", input.SLURM_job_name),
        string("#SBATCH --account=", input.SLURM_account_name),
        string("#SBATCH --partition=", input.SLURM_partition_name),
        string("#SBATCH --nodes=", input.SLURM_nodes_per_array_job),
        string("#SBATCH --ntasks=", input.SLURM_tasks_per_node),
        string("#SBATCH --cpus-per-task=", input.SLURM_cpus_per_task),
        string("#SBATCH --time=", input.SLURM_time_limit_dd_hhmmss),
        string("module load ", SLURM_module_load_R_version_name),
        string("time julia ", joinpath(run_directory_name, "run.jl \$SLURM_ARRAY_TASK_ID")),
    ]
    open(joinpath(run_directory_name, "run.slurm"), "w") do file
        write(file, join(slurm_script, "\n") * "\n")
    end
    # Submit the array of jobs
    # x = "-n"; y = "15"; z = "README.md"
    # P = Cmd([string("head"), x, y, z])
    # Q = Cmd(["tail", "-n", "5"])
    # run(P)
    # run(pipeline(P, Q))
    SHELL_COMMAND = Cmd([
        "sbatch",
        string("--array=1-", n_array_jobs, "%", input.SLURM_max_array_jobs_running),
        joinpath(run_directory_name, "run.jl"),
    ])
    run(SHELL_COMMAND)
    # Output
    nothing
end
