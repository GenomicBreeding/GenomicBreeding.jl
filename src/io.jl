"""
    mutable struct GBInput
        fname_geno::String
        fname_pheno::String
        bulk_cv::Bool
        populations::Union{Nothing, String, Vector{String}}
        traits::Union{Nothing, String, Vector{String}}
        models::Any
        n_folds::Int64
        n_replications::Int64
        maf::Float64
        mtv::Float64
        n_iter::Int64,
        n_burnin::Int64,
        fname_out_prefix::String,
        verbose::Bool
    end

Input struct (belongs to GBCore.AbstractGB type)

- `fname_geno`: genotype file (see file format guide: **TODO:** {URL HERE})
- `fname_pheno`: phenotype file (see file format guide: **TODO:** {URL HERE})
- `bulk_cv`: perform cross-validation across all populations, i.e. disregard population grouping (Default = false)
- `populations`: include only these populations (Default = nothing which means include all populations)
- `traits`: include only these traits (Default = nothing which means include all traits)
- `models`: include these model functions (Default = [ridge, bayesa]; see models list: **TODO:** {URL HERE})
- `n_folds`: number of k partitions for k-fold cross-validation (Default = 5)
- `n_replications`: number of replications for repeated k-fold cross-validation (Default = 5)
- `keep_all`: keep all entries upon merging genomes and phenomes potentially resulting in sparsities in both structs? (Default = false)
- `maf`: minimum allele frequency (Default = 0.05)
- `mtv`: minimum trait variance (Default = 1e-7)
- `n_iter`: number of Bayesian model fitting MCMC/HMC iteration (Default = 1_500)
- `n_burnin`: number of initial Bayesian model fitting MCMC/HMC iterations to be excluded from the posterior distribution (Default = 500)
- `fname_out_prefix`: prefix of the output files which may include directory names (Default = "" which translates to `./GBOuput/output-<yyyymmddHHMMSS>-<3_digit_random_number>-`)
- `SLURM_job_name`: name of the Slurm job array (Default = "")
- `SLURM_account_name`: Slurm account name (Default = "")
- `SLURM_partition_name`: Slurm paritition to use (Default = "")
- `SLURM_nodes_per_array_job`: number of nodes per array job (Default = 1)
- `SLURM_tasks_per_node`: number of tasks per node (Default = 1)
- `SLURM_cpus_per_task`: number of CPU cores per task (Default = 1)
- `SLURM_mem_G`: maximum memroy requested in gigabytes (Default = 1)
- `SLURM_time_limit_dd_hhmmss`: maximum computation time requested in the follwowing format: "dd-hh:mm:ss" (Default = "00-01:00:00")
- `SLURM_max_array_jobs_running`: maximum number of array jobs which can run simultaneously (Default = 20)
- `SLURM_module_load_R_version_name`: name of the R statistical language module name and version to be used to call R::BGLR (Default = "R")
- `verbose`: show messages (Default = true)
"""
mutable struct GBInput <: AbstractGB
    fname_geno::String
    fname_pheno::String
    bulk_cv::Bool
    populations::Union{Nothing,Vector{String}}
    traits::Union{Nothing,Vector{String}}
    models::Any
    n_folds::Int64
    n_replications::Int64
    keep_all::Bool
    maf::Float64
    mtv::Float64
    n_iter::Int64
    n_burnin::Int64
    fname_out_prefix::String
    SLURM_job_name::String
    SLURM_account_name::String
    SLURM_partition_name::String
    SLURM_nodes_per_array_job::Int64
    SLURM_tasks_per_node::Int64
    SLURM_cpus_per_task::Int64
    SLURM_mem_G::Int64
    SLURM_time_limit_dd_hhmmss::String
    SLURM_max_array_jobs_running::Int64
    SLURM_module_load_R_version_name::String
    verbose::Bool
    function GBInput(;
        fname_geno::String,
        fname_pheno::String,
        bulk_cv::Bool = false,
        populations::Union{Nothing,Vector{String}} = nothing,
        traits::Union{Nothing,Vector{String}} = nothing,
        models::Any = [ridge, bayesa],
        n_folds::Int64 = 5,
        n_replications::Int64 = 5,
        keep_all::Bool = false,
        maf::Float64 = 0.05,
        mtv::Float64 = 1e-7,
        n_iter::Int64 = 1_500,
        n_burnin::Int64 = 500,
        fname_out_prefix::String = "",
        SLURM_job_name::String = "",
        SLURM_account_name::String = "",
        SLURM_partition_name::String = "",
        SLURM_nodes_per_array_job::Int64 = 1,
        SLURM_tasks_per_node::Int64 = 1,
        SLURM_cpus_per_task::Int64 = 1,
        SLURM_mem_G::Int64 = 1,
        SLURM_time_limit_dd_hhmmss::String = "00-01:00:00",
        SLURM_max_array_jobs_running::Int64 = 20,
        SLURM_module_load_R_version_name::String = "R",
        verbose::Bool = true,
    )
        date = Dates.format(now(), "yyyymmddHHMMSS")
        randnum = Int(round(1_000 * rand()))
        fname_out_prefix = if fname_out_prefix == ""
            string("GBOutput/output-", date, "-", randnum, "-")
        else
            fname_out_prefix
        end
        SLURM_job_name = if SLURM_job_name == ""
            string("GBJob-", date, "-", randnum)
        else
            SLURM_job_name
        end
        new(
            fname_geno,
            fname_pheno,
            bulk_cv,
            populations,
            traits,
            models,
            n_folds,
            n_replications,
            keep_all,
            maf,
            mtv,
            n_iter,
            n_burnin,
            fname_out_prefix,
            SLURM_job_name,
            SLURM_account_name,
            SLURM_partition_name,
            SLURM_nodes_per_array_job,
            SLURM_tasks_per_node,
            SLURM_cpus_per_task,
            SLURM_mem_G,
            SLURM_time_limit_dd_hhmmss,
            SLURM_max_array_jobs_running,
            SLURM_module_load_R_version_name,
            verbose,
        )
    end
end

"""
    Base.hash(x::GBInput, h::UInt)::UInt

Hash a GBInput struct using all its fields.
We deliberately excluded the allele_frequencies, and mask for efficiency.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> input = GBInput(fname_geno="", fname_pheno="");

julia> typeof(hash(input))
UInt64
```
"""
function Base.hash(x::GBInput, h::UInt)::UInt
    # x = GBInput(fname_geno="", fname_pheno=""); h = 1.00
    hash(
        x.fname_geno,
        hash(
            x.fname_pheno,
            hash(
                x.bulk_cv,
                hash(
                    x.populations,
                    hash(
                        x.traits,
                        hash(
                            x.models,
                            hash(
                                x.n_folds,
                                hash(
                                    x.n_replications,
                                    hash(
                                        x.keep_all,
                                        hash(
                                            x.maf,
                                            hash(
                                                x.mtv,
                                                hash(
                                                    x.n_iter,
                                                    hash(
                                                        x.n_burnin,
                                                        hash(
                                                            x.fname_out_prefix,
                                                            hash(
                                                                x.SLURM_job_name,
                                                                hash(
                                                                    x.SLURM_account_name,
                                                                    hash(
                                                                        x.SLURM_partition_name,
                                                                        hash(
                                                                            x.SLURM_nodes_per_array_job,
                                                                            hash(
                                                                                x.SLURM_tasks_per_node,
                                                                                hash(
                                                                                    x.SLURM_cpus_per_task,
                                                                                    hash(
                                                                                        x.SLURM_mem_G,
                                                                                        hash(
                                                                                            x.SLURM_time_limit_dd_hhmmss,
                                                                                            hash(
                                                                                                x.SLURM_max_array_jobs_running,
                                                                                                hash(
                                                                                                    x.SLURM_module_load_R_version_name,
                                                                                                    hash(x.verbose, h),
                                                                                                ),
                                                                                            ),
                                                                                        ),
                                                                                    ),
                                                                                ),
                                                                            ),
                                                                        ),
                                                                    ),
                                                                ),
                                                            ),
                                                        ),
                                                    ),
                                                ),
                                            ),
                                        ),
                                    ),
                                ),
                            ),
                        ),
                    ),
                ),
            ),
        ),
    )
end


"""
    Base.:(==)(x::GBInput, y::GBInput)::Bool

Equality of GBInput structs using the hash function defined for GBInput structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> input_1 = input = GBInput(fname_geno="geno1.jld2", fname_pheno="pheno1.jld2", fname_out_prefix="test1-", SLURM_job_name="slurmjob1");

julia> input_2 = input = GBInput(fname_geno="geno1.jld2", fname_pheno="pheno1.jld2", fname_out_prefix="test1-", SLURM_job_name="slurmjob1");

julia> input_3 = input = GBInput(fname_geno="geno2.jld2", fname_pheno="pheno2.jld2");

julia> input_1 == input_2
true

julia> input_1 == input_3
false
```
"""
function Base.:(==)(x::GBInput, y::GBInput)::Bool
    hash(x) == hash(y)
end

"""
    clone(x::GBInput)::GBInput

Clone a GBInput object

## Example
```jldoctest; setup = :(using GenomicBreeding)
julia> input = GBInput(fname_geno="geno1.jld2", fname_pheno="pheno1.jld2");

julia> copy_input = clone(input);

julia> input == copy_input
true
```
"""
function clone(x::GBInput)::GBInput
    GBInput(
        fname_geno = deepcopy(x.fname_geno),
        fname_pheno = deepcopy(x.fname_pheno),
        bulk_cv = deepcopy(x.bulk_cv),
        populations = deepcopy(x.populations),
        traits = deepcopy(x.traits),
        models = deepcopy(x.models),
        n_folds = deepcopy(x.n_folds),
        n_replications = deepcopy(x.n_replications),
        keep_all = deepcopy(x.keep_all),
        maf = deepcopy(x.maf),
        mtv = deepcopy(x.mtv),
        n_iter = deepcopy(x.n_iter),
        n_burnin = deepcopy(x.n_burnin),
        fname_out_prefix = deepcopy(x.fname_out_prefix),
        SLURM_job_name = deepcopy(x.SLURM_job_name),
        SLURM_account_name = deepcopy(x.SLURM_account_name),
        SLURM_partition_name = deepcopy(x.SLURM_partition_name),
        SLURM_nodes_per_array_job = deepcopy(x.SLURM_nodes_per_array_job),
        SLURM_tasks_per_node = deepcopy(x.SLURM_tasks_per_node),
        SLURM_cpus_per_task = deepcopy(x.SLURM_cpus_per_task),
        SLURM_mem_G = deepcopy(x.SLURM_mem_G),
        SLURM_time_limit_dd_hhmmss = deepcopy(x.SLURM_time_limit_dd_hhmmss),
        SLURM_max_array_jobs_running = deepcopy(x.SLURM_max_array_jobs_running),
        SLURM_module_load_R_version_name = deepcopy(x.SLURM_module_load_R_version_name),
        verbose = deepcopy(x.verbose),
    )
end

"""
    checkdims(input::GBInput)::Bool

Check dimension compatibility of the fields of the GBInput struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> input = GBInput(fname_geno="geno1.jld2", fname_pheno="pheno1.jld2");

julia> checkdims(input)
true

julia> input.models = nothing

julia> checkdims(input)
false
```
"""
function GBCore.checkdims(input::GBInput)::Bool
    !isnothing(input.models)
end

"""
    load(input::GBInput)::Tuple{Genomes, Phenomes}

Load, merge and filter genotype and phenotype data

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    
julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], verbose=false);

julia> genomes, phenomes = load(input);

julia> length(unique(genomes.populations)) == length(unique(phenomes.populations)) == 2
true

julia> length(phenomes.traits) == 1
true

julia> rm.([fname_geno, fname_pheno]);
```
"""
function load(input::GBInput)::Tuple{Genomes,Phenomes}
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = writedelimited(genomes, fname="test-geno.tsv")
    # # fname_geno = writejld2(genomes, fname="test-geno.jld2")
    # # fname_geno = writevcf(genomes, fname="test-geno.vcf") ### Will have unknown population groupings
    # fname_pheno = writedelimited(phenomes, fname="test-pheno.tsv")
    # # fname_pheno = writejld2(phenomes, fname="test-pheno.jld2")
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # Parse input
    fname_geno = input.fname_geno
    fname_pheno = input.fname_pheno
    bulk_cv = input.bulk_cv
    populations = input.populations
    traits = input.traits
    n_folds = input.n_folds
    keep_all = input.keep_all
    maf = input.maf
    mtv = input.mtv
    verbose = input.verbose
    # Load genomes and phenomes
    genomes = try
        readdelimited(Genomes, fname = fname_geno)
    catch
        try
            readjld2(Genomes, fname = fname_geno)
        catch
            try
                readvcf(fname = fname_geno)
                @warn "The vcf file: `" *
                      fname_geno *
                      "` lacks population grouping information. This will likely result in conflicts in merging with th phenotype data."
            catch
                throw(
                    ArgumentError(
                        "Unrecognised genotype file format: `" *
                        fname_geno *
                        "`.\n" *
                        "Please refer to the file format guide: **TODO:** {URL HERE}",
                    ),
                )
            end
        end
    end
    phenomes = try
        readdelimited(Phenomes, fname = fname_pheno)
    catch
        try
            readjld2(Phenomes, fname = fname_pheno)
        catch
            throw(
                ArgumentError(
                    "Unrecognised phenotype file format: `" *
                    fname_pheno *
                    "`.\n" *
                    "Please refer to the file format guide: **TODO:** {URL HERE}",
                ),
            )
        end
    end
    if verbose
        println("##############")
        println("### Loaded ###")
        println("##############")
        println("- genomes:")
        display(dimensions(genomes))
        println("- phenomes:")
        display(dimensions(phenomes))
    end
    # Merge the genomes and phenomes
    genomes, phenomes = merge(genomes, phenomes, keep_all = keep_all)
    if verbose
        println("##############")
        println("### Merged ###")
        println("##############")
        println("- genomes:")
        display(dimensions(genomes))
        println("- phenomes:")
        display(dimensions(phenomes))
    end
    # Prepare population/s and trait/s to retain
    populations = if isnothing(populations)
        unique(sort(phenomes.populations))
    else
        populations
    end
    traits = if isnothing(traits)
        sort(phenomes.traits)
    else
        traits
    end
    # Check population names
    all_populations = unique(sort(phenomes.populations))
    unrecognised_populations = []
    for pop in populations
        if !(pop ∈ all_populations)
            push!(unrecognised_populations, pop)
        end
    end
    if length(unrecognised_populations) > 0
        throw(ArgumentError("Unrecognised population/s:\n\t‣ " * join(unrecognised_populations, "\n\t‣ ")))
    end
    # Check trait names
    unrecognised_traits = []
    for trait in traits
        if !(trait ∈ phenomes.traits)
            push!(unrecognised_traits, trait)
        end
    end
    if length(unrecognised_traits) > 0
        throw(ArgumentError("Unrecognised trait/s:\n\t‣ " * join(unrecognised_traits, "\n\t‣ ")))
    end
    # Define the entries to retain in both phenomes and genomes
    idx_entries = findall([sum(populations .== pop) > 0 for pop in genomes.populations])
    # Slice phenomes to retain only the requested population/s and/or trait/s
    if length(populations) < length(all_populations)
        phenomes = slice(phenomes, idx_entries = idx_entries)
    end
    if length(traits) < length(phenomes.traits)
        idx_traits = findall([sum(traits .== trait) > 0 for trait in phenomes.traits])
        phenomes = slice(phenomes, idx_traits = idx_traits)
    end
    # Filter phenomes by mtv, i.e. retain traits with variance
    fixed_traits = []
    for trait in phenomes.traits
        # trait = phenomes.traits[1]
        y = phenomes.phenotypes[:, phenomes.traits.==trait][:, 1]
        y = y[.!ismissing.(y).&&.!isnan.(y).&&.!isinf.(y)]
        if length(y) > 0
            if var(y) < mtv
                push!(fixed_traits, trait)
            end
        else
            push!(fixed_traits, trait)
        end
    end
    if length(fixed_traits) > 0
        @warn "Excluding the following fixed (or all values are missing/NaN/Infinities) trait/s:\n\t‣ " *
              join(fixed_traits, "\n\t‣ ")
        idx_traits = findall([sum(fixed_traits .== trait) == 0 for trait in phenomes.traits])
        if length(idx_traits) == 0
            throw(
                ArgumentError(
                    "All requested traits are fixed, i.e.\n\t‣ " *
                    join(phenomes.traits, "\n\t‣ ") *
                    "\n. On the following requested populations:\n\t‣ " *
                    join(populations, "\n\t‣ "),
                ),
            )
        end
        phenomes = slice(phenomes, idx_traits = idx_traits)
    end
    # Check if we have enough entries per fold
    if !bulk_cv && (length(unique(phenomes.populations)) > 1)
        # Per trait per population, is we do not wish bulk CV and we have more than 1 population
        for trait in phenomes.traits
            for population in phenomes.populations
                # trait=phenomes.traits[1]; population=phenomes.populations[1];
                y = phenomes.phenotypes[phenomes.populations.==population, phenomes.traits.==trait][:, 1]
                y = y[.!ismissing.(y).&&.!isnan.(y).&&.!isinf.(y)]
                if ceil(length(y) / n_folds) < 5
                    throw(
                        ArgumentError(
                            "The numer of entries in population: `" *
                            population *
                            "` (n=" *
                            string(length(y)) *
                            ") results in less than 5 entries per " *
                            string(n_folds) *
                            "-fold cross-validation. " *
                            "Please consider reducing the number of " *
                            "folds for this specific population or setting `cv_bulk = true`.",
                        ),
                    )
                end
            end
        end
    else
        # Per trait across all populations - bulked
        for trait in phenomes.traits
            # trait = phenomes.traits[1]
            y = phenomes.phenotypes[:, phenomes.traits.==trait][:, 1]
            y = y[.!ismissing.(y).&&.!isnan.(y).&&.!isinf.(y)]
            if ceil(length(y) / n_folds) < 5
                throw(
                    ArgumentError(
                        "The numer of entries across all populations " *
                        "results in less than 5 entries per " *
                        string(n_folds) *
                        "-fold cross-validation. " *
                        "Please consider reducing the number of folds or " *
                        "setting `cv_bulk = true`.",
                    ),
                )
            end
        end
    end
    # Finally, slice the genomes to retain only the requested populations
    # Note that we delayed this for computational efficiency,
    # i.e. to catch errors with the Phenomes struct before processing this potentially massive Genomes struct.
    if length(populations) < length(all_populations)
        genomes = slice(genomes, idx_entries = idx_entries)
    end
    # Filter genomes by maf (minimum allele frequency)
    genomes = filter(genomes, maf)
    # Show dimensions of the input genomes and phenomes after merging and filterings
    if verbose
        println("################")
        println("### Filtered ###")
        println("################")
        println("- genomes:")
        display(dimensions(genomes))
        println("- phenomes:")
        display(dimensions(phenomes))
    end
    # Output
    if genomes.entries != phenomes.entries
        throw(
            ErrorException(
                "Error loading the genomes and phenomes data from:\n\t‣ Genotype file: " *
                fname_geno *
                "\n\t‣ Phenotype file: " *
                fname_pheno,
            ),
        )
    end
    (genomes, phenomes)
end


"""
    prepareinputs(input::GBInput)::Vector{GBInput}

Prepare GBInputs for Slurm array jobs

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    
julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, verbose=false);

julia> inputs = prepareinputs(input);

julia> length(inputs) == 30
true

julia> rm.([fname_geno, fname_pheno]);
```
"""
function prepareinputs(input::GBInput)::Vector{GBInput}
    # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # phenomes = extractphenomes(trials)
    # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno)
    # Load genomes and phenomes to check their validity and dimensions
    genomes, phenomes = load(input)
    # Count the number of models, traits, and populations
    m = length(input.models)
    t = length(phenomes.traits)
    p = length(unique(phenomes.populations))
    # Prepare the GBInputs
    inputs = if p > 1
        Vector{GBInput}(undef, m * t * (p + 2))
    else
        Vector{GBInput}(undef, m * t * p)
    end
    i = 1
    for model in input.models
        # model = input.models[1]
        for trait in phenomes.traits
            # trait = phenomes.traits[1]
            populations = if p > 1
                vcat("BULK_CV", "ACROSS_POP_CV", sort(unique(phenomes.populations)))
            else
                sort(unique(phenomes.populations))
            end
            for population in populations
                # population = populations[1]
                bulk_cv, population = if population == "BULK_CV"
                    true, nothing
                elseif population == "ACROSS_POP_CV"
                    false, populations[3:end]
                else
                    false, [population]
                end
                inputs[i] = GBInput(
                    fname_geno = input.fname_geno,
                    fname_pheno = input.fname_pheno,
                    bulk_cv = bulk_cv,
                    populations = population,
                    traits = [trait],
                    models = [model],
                    n_folds = input.n_folds,
                    n_replications = input.n_replications,
                    keep_all = input.keep_all,
                    maf = input.maf,
                    mtv = input.mtv,
                    n_iter = input.n_iter,
                    n_burnin = input.n_burnin,
                    fname_out_prefix = input.fname_out_prefix,
                    verbose = input.verbose,
                )
                i += 1
            end
        end
    end
    # Output
    inputs
end
