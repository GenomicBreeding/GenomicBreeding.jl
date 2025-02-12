"""
    submitslurmarrayjobs(; input::GBInput, analysis::Function)::String

Assess genomic prediction accuracy via replicated k-fold cross-validation.
Outputs are saved as JLD2 (each containing a CV struct per fold, replication, and trait) and possibly text file/s containing notes describing why some jobs failed.

# Example
<!-- ```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames) -->
```
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5, fname_out_prefix="GBOutput/test-", verbose=false);

julia> outdir = submitslurmarrayjobs(input=input, analysis=assess)
GBOutput

julia> run(`squeue`)
```
"""
function submitslurmarrayjobs(; input::GBInput, analysis::Function)::String
    # genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=9)
    # analysis = assess
    # Catch the 3 main errors early:
    #   1. input genotype file does not exist
    #   2. input phenotype file does not exist
    #   3. R::BGLR package is not installed.
    errors::Vector{String} = []
    if !isfile(input.fname_geno)
        push!(errors, string("The genotype file (", input.fname_geno, ") does not exist."))
    end
    if !isfile(input.fname_pheno)
        push!(errors, string("The phenotype file (", input.fname_pheno, ") does not exist."))
    end
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
        writejld2(inputs[i], fname = joinpath(run_outdir, string("GBInput-", i, ".jld2")))
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
        "import GenomicBreeding: ols, ridge, lasso, bayesa, bayesb, bayesc",
        string("input = readjld2(GBInput, fname=joinpath(\"", run_outdir, "\", string(\"GBInput-\", i, \".jld2\")))"),
        "display(input)",
        string("output = ", analysis, "(input)"),
        "display(output)",
    ]
    open(joinpath(run_outdir, "run.jl"), "w") do file
        write(file, join(julia_script, "\n") * "\n")
    end
    # Save the Slurm run file
    slurm_script = [
        "#!/bin/bash",
        string("#SBATCH --job-name='", input.SLURM_job_name, "'"),
        string("#SBATCH --account='", input.SLURM_account_name, "'"),
        string("#SBATCH --partition='", input.SLURM_partition_name, "'"),
        string("#SBATCH --nodes=", input.SLURM_nodes_per_array_job),
        string("#SBATCH --ntasks=", input.SLURM_tasks_per_node),
        string("#SBATCH --cpus-per-task=", input.SLURM_cpus_per_task),
        string("#SBATCH --mem=", input.SLURM_mem_G, "G"),
        string("#SBATCH --time=", input.SLURM_time_limit_dd_hhmmss),
        string("module load ", input.SLURM_module_load_R_version_name),
        "LD_LIBRARY_PATH=\"\"",
        string(
            "time julia --threads ",
            input.SLURM_cpus_per_task,
            " ",
            joinpath(run_outdir, "run.jl \$SLURM_ARRAY_TASK_ID"),
        ),
    ]
    if input.SLURM_account_name == ""
        slurm_script = slurm_script[isnothing.(match.(Regex("^#SBATCH --account="), slurm_script))]
    end
    if input.SLURM_partition_name == ""
        slurm_script = slurm_script[isnothing.(match.(Regex("^#SBATCH --partition="), slurm_script))]
    end
    open(joinpath(run_outdir, "run.slurm"), "w") do file
        write(file, join(slurm_script, "\n") * "\n")
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

# using Libdl
# filter!(contains("curl"), dllist())