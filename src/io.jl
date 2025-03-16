"""
    mutable struct GBInput

Main input struct for genomic breeding analysis (implements GBCore.AbstractGB)

# Fields

## Required
- `fname_geno`: Path to genotype file (see file format guide for details)

## Optional Data Files
- `fname_pheno`: Path to phenotype file (Default: "" - see file format guide for details)
- `fname_allele_effects_jld2s`: Vector of paths to JLD2 files containing Fit structs (Default: [""])

## Analysis Settings
- `analysis`: Analysis function to perform (Default: cv)
  - `cv`: Replicated k-fold cross-validation
  - `fit`: Fit genomic prediction models to extract allele effects
  - `predict`: Compute GEBVs using existing model fits
  - `gwas`: Genome-wide association study
- `bulk_cv`: Perform cross-validation across all populations (Default: false)
- `populations`: Vector of populations to include (Default: all)
- `traits`: Vector of traits to analyze (Default: all)
- `models`: Vector of genomic prediction model functions (Default: [ridge, bayesa])
- `gwas_models`: Vector of GWAS model functions (Default: [gwasols, gwaslmm])

## Cross-Validation Parameters
- `n_folds`: Number of CV folds (Default: 5)
- `n_replications`: Number of CV replications (Default: 5)

## Filtering Parameters
- `keep_all`: Keep all entries when merging data (Default: false)
- `maf`: Minimum allele frequency (Default: 0.05)
- `mtv`: Minimum trait variance (Default: 1e-7)

## Model Parameters
- `n_iter`: MCMC iterations (Default: 1_500)
- `n_burnin`: MCMC burn-in iterations (Default: 500)

## Output Settings
- `fname_out_prefix`: Output file prefix & directory (Default: GBOutput/output-<timestamp>-)
- `verbose`: Show progress messages (Default: true)

## SLURM Settings
- `SLURM_job_name`: Job array name (Default: GBJob-<timestamp>)
- `SLURM_account_name`: Account name (Default: "")
- `SLURM_partition_name`: Partition to use (Default: "")
- `SLURM_nodes_per_array_job`: Nodes per job (Default: 1)
- `SLURM_tasks_per_node`: Tasks per node (Default: 1)
- `SLURM_cpus_per_task`: CPUs per task (Default: 1)
- `SLURM_mem_G`: Memory in GB (Default: 1)
- `SLURM_time_limit_dd_hhmmss`: Time limit as "dd-hh:mm:ss" (Default: "00-01:00:00")
- `SLURM_max_array_jobs_running`: Max concurrent array jobs (Default: 20)
- `SLURM_module_load_Conda_version_name`: Conda module name (Default: "Miniconda3")
- `SLURM_module_load_R_version_name`: R module name (Default: "conda" which will use the R installed in the conda environment - see installation instructions for details)
- `SLURM_module_load_Julia_version_name`: Julia module name (Default: "" which will use the Julia installed via JuliaUp - see installation instructions for details)
"""
mutable struct GBInput <: AbstractGB
    fname_geno::String
    fname_pheno::String
    fname_allele_effects_jld2s::Vector{String}
    analysis::Function
    bulk_cv::Bool
    populations::Union{Nothing,Vector{String}}
    traits::Union{Nothing,Vector{String}}
    models::Any
    n_folds::Int64
    n_replications::Int64
    gwas_models::Any
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
    SLURM_module_load_Conda_version_name::String
    SLURM_module_load_R_version_name::String
    SLURM_module_load_Julia_version_name::String
    verbose::Bool
    function GBInput(;
        fname_geno::String,
        fname_pheno::String = "",
        fname_allele_effects_jld2s::Vector{String} = [""],
        analysis::Function = cv,
        bulk_cv::Bool = false,
        populations::Union{Nothing,Vector{String}} = nothing,
        traits::Union{Nothing,Vector{String}} = nothing,
        models::Any = [ridge, bayesa],
        n_folds::Int64 = 5,
        n_replications::Int64 = 5,
        gwas_models::Any = [gwasols, gwaslmm],
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
        SLURM_module_load_Conda_version_name::String = "Miniconda3",
        SLURM_module_load_R_version_name::String = "conda",
        SLURM_module_load_Julia_version_name::String = "",
        verbose::Bool = true,
    )
        date = Dates.format(now(), "yyyymmddHHMMSSssss")
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
            fname_allele_effects_jld2s,
            analysis,
            bulk_cv,
            populations,
            traits,
            models,
            n_folds,
            n_replications,
            gwas_models,
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
            Float64(SLURM_mem_G),
            SLURM_time_limit_dd_hhmmss,
            SLURM_max_array_jobs_running,
            SLURM_module_load_Conda_version_name,
            SLURM_module_load_R_version_name,
            SLURM_module_load_Julia_version_name,
            verbose,
        )
    end
end

"""
    Base.hash(x::GBInput, h::UInt)::UInt

Compute a hash value for a `GBInput` struct by combining the hash values of all its fields.

# Arguments
- `x::GBInput`: The input structure to be hashed
- `h::UInt`: The hash value seed

# Returns
- `UInt`: A hash value that uniquely identifies the content of the `GBInput` struct

# Details
This method implements hash computation for the `GBInput` type by iterating through
all fields and combining their hash values.

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> input = GBInput(fname_geno="", fname_pheno="");

julia> typeof(hash(input))
UInt64
```
"""
function Base.hash(x::GBInput, h::UInt)::UInt
    # x = GBInput(fname_geno="", fname_pheno=""); h = 1.00
    for field in fieldnames(typeof(x))
        # field = fieldnames(typeof(x))[1]
        h = hash(getfield(x, field), h)
    end
    h
end


"""
    Base.:(==)(x::GBInput, y::GBInput)::Bool

Compare two `GBInput` structs for equality by comparing their hash values.

This method overloads the `==` operator for `GBInput` structs, allowing direct comparison
using the `==` operator. Two `GBInput` structs are considered equal if they have identical
hash values, which implies they have the same values for all relevant fields.

# Arguments
- `x::GBInput`: First GBInput struct to compare
- `y::GBInput`: Second GBInput struct to compare

# Returns
- `Bool`: `true` if the hash values of both structs are equal, `false` otherwise

# Examples
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

Create a deep copy of a GBInput object, duplicating all field values into a new instance.

This function allows you to create an independent copy of a GBInput object where modifications 
to the clone won't affect the original object.

# Arguments
- `x::GBInput`: The source GBInput object to be cloned

# Returns
- `::GBInput`: A new GBInput instance with identical field values

# Example
```jldoctest; setup = :(using GenomicBreeding)
julia> input = GBInput(fname_geno="geno1.jld2", fname_pheno="pheno1.jld2");

julia> copy_input = clone(input);

julia> input == copy_input
true
```
"""
function GBCore.clone(x::GBInput)::GBInput
    # x = GBInput(fname_geno=""); x.fname_geno="test_geno.jld2"; x.fname_pheno = "test_pheno.tsv";
    out = GBInput(fname_geno = "", fname_pheno = "")
    for field in fieldnames(typeof(x))
        # field = fieldnames(typeof(x))[1]
        setfield!(out, field, getfield(x, field))
    end
    out
end

"""
    checkdims(input::GBInput)::Bool

Check dimension compatibility of the fields of the GBInput struct.
Returns `true` if both `models` and `gwas_models` fields in the input are not `nothing`,
`false` otherwise.

# Arguments
- `input::GBInput`: Input structure containing genomic data and models

# Returns
- `Bool`: `true` if dimensions are compatible, `false` otherwise

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
    !isnothing(input.models) & !isnothing(input.gwas_models)
end

"""
    checkinputs(input::GBInput)::Vector{String}

Check the compatibility and validity of inputs for genomic analysis.

Returns a vector of error messages. An empty vector indicates all inputs are valid.

# Arguments
- `input::GBInput`: A struct containing analysis parameters including:
  - `analysis`: Analysis type (`cv`, `fit`, `predict`, or `gwas`)
  - `fname_geno`: Path to genotype file
  - `fname_pheno`: Path to phenotype file
  - `models`: Vector of selected models
  - `fname_allele_effects_jld2s`: Paths to allele effects files (for predict)

# Returns
- `Vector{String}`: Collection of error messages, empty if all inputs are valid

# Validation Rules
- Analysis type must be one of: `cv`, `fit`, `predict`, or `gwas`
- For `cv`, `fit`, and `gwas`:
  - Genotype file must exist
  - Phenotype file must exist
  - At least one model must be specified
- For `predict`:
  - Genotype file must exist
  - All specified allele effects files must exist

# Examples
```jldoctest; setup = :(using GBCore, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    
julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=cv);

julia> length(checkinputs(input)) == 0
true

julia> input.fname_pheno = ""; length(checkinputs(input)) == 0
false
```
"""
function checkinputs(input::GBInput)::Vector{String}
    errors::Vector{String} = []
    valid_analysis_functions = [cv, fit, predict, gwas]
    if !(input.analysis ∈ valid_analysis_functions)
        push!(
            errors,
            string(
                "Analysis: `",
                string(input.analysis),
                "` invalid. Please choose from:\n\t‣ ",
                join(string.(valid_analysis_functions), "\n\t‣ "),
            ),
        )
    end
    if input.analysis ∈ [cv, fit, gwas]
        if !isfile(input.fname_geno)
            push!(errors, string("The input genotype file `", input.fname_geno, "` does not exist."))
        end
        if !isfile(input.fname_pheno)
            push!(errors, string("The input phenotype file `", input.fname_pheno, "` does not exist."))
        end
        if length(input.models) == 0
            if input.analysis ∈ [cv, fit]
                push!(
                    errors,
                    string(
                        "No model specified for `",
                        input.analysis,
                        "`. Please specify the genomic prediction model/s you wish to use, e.g. `rigde`, `bayesa` or `bayesb`.",
                    ),
                )
            else
                push!(
                    errors,
                    string(
                        "No model specified for `",
                        input.analysis,
                        "`. Please specify the genomic prediction model/s you wish to use, e.g. `gwasols`, `gwaslmm` or `gwasreml`.",
                    ),
                )
            end
        end
    end
    if input.analysis ∈ [predict]
        if !isfile(input.fname_geno)
            push!(errors, string("The input genotype file `", input.fname_geno, "` does not exist."))
        end
        for fname in input.fname_allele_effects_jld2s
            if !isfile(fname)
                push!(
                    errors,
                    string("The input allele effects file `", input.fname_allele_effects_jld2s, "` does not exist."),
                )
            end
        end
    end
    return errors
end


"""
    loadgenomesphenomes(input::GBInput)::Tuple{Genomes, Phenomes, Vector{String}, Vector{String}}

Load, merge, and filter genotype and phenotype data from specified input files.

# Arguments
- `input::GBInput`: A struct containing input parameters including:
  - `fname_geno`: Path to genotype data file
  - `fname_pheno`: Path to phenotype data file
  - `bulk_cv`: Boolean for bulk cross-validation
  - `populations`: Vector of population names to include (optional)
  - `traits`: Vector of trait names to include (optional)
  - `n_folds`: Number of cross-validation folds
  - `keep_all`: Boolean to keep all entries during merging
  - `maf`: Minimum allele frequency threshold
  - `mtv`: Minimum trait variance threshold
  - `verbose`: Boolean for detailed output

# Returns
A tuple containing:
1. `Genomes`: Filtered genomic data
2. `Phenomes`: Filtered phenotypic data
3. `Vector{String}`: Traits to skip due to insufficient data for cross-validation
4. `Vector{String}`: Populations to skip due to insufficient data for cross-validation

# Notes
- Supports multiple file formats (string-delimited, JLD2, VCF)
- Performs data validation and compatibility checks
- Filters traits with variance below minimum trait variance (mtv) threshold
- Ensures sufficient sample size for cross-validation
- Filters markers based on minimum allele frequency (maf)

# Examples
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    
julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, populations=["pop_1", "pop_3"], traits=["trait_1"], verbose=false);

julia> genomes, phenomes, traits_to_skip, populations_to_skip = loadgenomesphenomes(input);

julia> length(unique(genomes.populations)) == length(unique(phenomes.populations)) == 2
true

julia> length(phenomes.traits) == 1
true

julia> rm.([fname_geno, fname_pheno]);
```
"""
function loadgenomesphenomes(input::GBInput)::Tuple{Genomes,Phenomes,Vector{String},Vector{String}}
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
    # Check input files and analysis compatibility
    errors = checkinputs(input)
    if length(errors) > 0
        throw(ArgumentError(string("Errors:\n\t• ", join(errors, "\n\t• "))))
    end
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
    phenomes = if (fname_pheno != "") && (input.analysis ∈ [cv, fit, gwas])
        try
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
                        "You may have loaded a Trials file. " *
                        "Please refer to the file format guide: **TODO:** {URL HERE}",
                    ),
                )
            end
        end
    else
        dummpy_phenomes = Phenomes(n = length(genomes.entries), t = 1)
        dummpy_phenomes.entries = genomes.entries
        dummpy_phenomes.populations = genomes.populations
        dummpy_phenomes.traits = ["dummy_trait"]
        dummpy_phenomes.phenotypes = rand(length(dummpy_phenomes.entries), 1)
        dummpy_phenomes
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
    traits = if isnothing(traits) || (traits == [""])
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
    # Trait-population combinations to skip
    traits_to_skip::Vector{String} = []
    populations_to_skip::Vector{String} = []
    # Check if we have enough entries per fold
    if !bulk_cv && (length(unique(phenomes.populations)) > 1)
        # Per trait per population, is we do not wish bulk CV and we have more than 1 population
        for trait in phenomes.traits
            for population in phenomes.populations
                # trait=phenomes.traits[16]; population=phenomes.populations[16];
                y = phenomes.phenotypes[phenomes.populations.==population, phenomes.traits.==trait][:, 1]
                y = y[.!ismissing.(y).&&.!isnan.(y).&&.!isinf.(y)]
                if ceil(length(y) / n_folds) < 5
                    # throw(
                    #     ArgumentError(
                    #         "The numer of entries in population: `" *
                    #         population *
                    #         "` (n=" *
                    #         string(length(y)) *
                    #         ") results in less than 5 entries per " *
                    #         string(n_folds) *
                    #         "-fold cross-validation. " *
                    #         "Please consider reducing the number of " *
                    #         "folds for this specific population or setting `cv_bulk = true`.",
                    #     ),
                    # )
                    # println(string("trait = \"", trait, "\"; population = \"", population, "\"; ceil(length(y) / n_folds) = ", ceil(length(y) / n_folds)))
                    push!(traits_to_skip, trait)
                    push!(populations_to_skip, population)
                end


            end
        end
    else
        # Per trait across all populations - BULK_CV
        for trait in phenomes.traits
            # trait = phenomes.traits[16]
            y = phenomes.phenotypes[:, phenomes.traits.==trait][:, 1]
            y = y[.!ismissing.(y).&&.!isnan.(y).&&.!isinf.(y)]
            if ceil(length(y) / n_folds) < 5
                # throw(
                #     ArgumentError(
                #         "The numer of entries across all populations " *
                #         "results in less than 5 entries per " *
                #         string(n_folds) *
                #         "-fold cross-validation. " *
                #         "Please consider reducing the number of folds or " *
                #         "setting `cv_bulk = true`.",
                #     ),
                # )
                push!(traits_to_skip, trait)
                push!(populations_to_skip, "BULK_CV")
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
    (genomes, phenomes, traits_to_skip, populations_to_skip)
end

"""
    loadcvs(input::GBInput; min_train_size::Int64=10)::Vector{CV}

Load and filter cross-validation (CV) results from files generated by `GenomicBreeding.cv()`.

# Arguments
- `input::GBInput`: Input configuration containing the output directory path in `fname_out_prefix`
- `min_train_size::Int64=10`: Minimum required size for training sets (default: 10)

# Returns
- `Vector{CV}`: Vector of valid CV objects that meet the following criteria:
    - Contains complete metrics (length of fit.metrics == 9)
    - Has valid correlation metric (not missing, NaN, or Inf)
    - Validation set size ≥ min_train_size

# Throws
- `ArgumentError`: If the output directory doesn't exist or contains no CV results

# Details
The function searches for files with pattern "-cv-" and extension ".jld2" in the output directory.
Invalid or failed CV results are automatically filtered out during loading.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, fname_out_prefix="GBOutput_cv_1/output-1-", populations=["pop_1", "pop_3"], traits=["trait_1"], n_replications=2, n_folds=3, verbose=false);

julia> fnames_cvs, fnames_notes = cv(input);

julia> cvs = loadcvs(input);

julia> length(cvs) == length(fnames_cvs)
true
```
"""
function loadcvs(input::GBInput; min_train_size::Int64 = 10)::Vector{CV}
    directory_name = if dirname(input.fname_out_prefix) == ""
        pwd()
    else
        dirname(input.fname_out_prefix)
    end
    if !isdir(directory_name)
        throw(
            ArgumentError(
                "The output directory (`" *
                directory_name *
                "`) which should contain the output of `GenomicBreeding.cv(...)` does not exist. " *
                "Have you run `submitslurmarrayjobs(...)` using an input with `analysis=cv`?",
            ),
        )
    end
    files = readdir(directory_name)
    idx = findall(.!isnothing.(match.(Regex("-cv-"), files)) .&& .!isnothing.(match.(Regex("jld2\$"), files)))
    if length(idx) == 0
        throw(
            ArgumentError(
                "The output directory (`" *
                directory_name *
                "`) does not contain the output of `GenomicBreeding.cv(...)`. " *
                "Have you run `submitslurmarrayjobs(...)` using an input with `analysis=cv`?",
            ),
        )
    end
    fnames_cvs = joinpath.(directory_name, files[idx])
    # println(idx) 
    # println(fnames_cvs) 
    cvs::Vector{CV} = []
    for (i, fname) in enumerate(fnames_cvs)
        # i = 1; fname = fnames_cvs[i];
        try
            c = readjld2(CV, fname = fname)
            # Skip folds which failed due to insufficient sizes and/or nil phenotype variance and/or training size less than 10

            if (
                (length(c.fit.metrics) == 9) &&
                !ismissing(c.metrics["cor"]) &&
                !isnan(c.metrics["cor"]) &&
                !isinf(c.metrics["cor"]) &&
                (length(c.validation_entries) >= min_train_size)
            )
                push!(cvs, c)
            end
        catch
            # println(fname) 
            continue
        end
    end
    cvs
end

"""
    loadfits(input::GBInput)::Vector{Fit}

Load fitted allele frequency effects from genomic prediction models stored in JLD2 files.

# Arguments
- `input::GBInput`: Input configuration containing paths to JLD2 files with fitted allele effects.

# Returns
- `Vector{Fit}`: Array of `Fit` objects containing the loaded allele frequency effects.

# Details
The function attempts to load all JLD2 files specified in `input.fname_allele_effects_jld2s`.
If a file cannot be loaded, that entry will be skipped and remain `undef` in the output vector.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> proportion_of_variance = fill(0.0, 9, 3); proportion_of_variance[1, :] .= 1.00; # 100% variance on the additive genetic effects

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, proportion_of_variance=proportion_of_variance, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, fname_out_prefix="GBOutput_fit_2/output-2-", populations=["pop_1", "pop_3"], traits=["trait_1"], models=[bayesa, bayesb], verbose=false);

julia> input.fname_allele_effects_jld2s = GenomicBreeding.fit(input);

julia> fits = loadfits(input);

julia> length(fits) == length(input.fname_allele_effects_jld2s)
true
```
"""
function loadfits(input::GBInput)::Vector{Fit}
    fits = Vector{Fit}(undef, length(input.fname_allele_effects_jld2s))
    for (i, fname) in enumerate(input.fname_allele_effects_jld2s)
        # i = 1; fname = fname_allele_effects_jld2s[i];
        try
            fits[i] = readjld2(Fit, fname = fname)
        catch
            continue
        end
    end
    fits
end

"""
    prepareinputs(input::GBInput)::Vector{GBInput}

Create a vector of `GBInput` objects for parallel processing based on the input configuration.

# Arguments
- `input::GBInput`: Initial input configuration containing analysis parameters.

# Returns
- `Vector{GBInput}`: Array of `GBInput` objects configured for different combinations of:
  - Models (RR-BLUP, BayesB, etc.)
  - Traits from phenotype data
  - Population groups (including bulk and across-population analyses)

# Details
For different analysis types, the function generates the following combinations:
- Cross-validation (`cv`): 2 models × traits × (populations + bulk + across-pop)
- Model fitting (`fit`): 2 models × traits × (populations + bulk)
- Prediction (`predict`): 1 GBInput per allele effects file
- GWAS (`gwas`): 1 model × traits × (populations + bulk)

The function handles special cases:
- Skips trait-population combinations with insufficient data
- Adjusts settings for bulk and across-population analyses
- Configures specific parameters for prediction tasks

# Examples
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    
julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;

julia> input_cv = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=GenomicBreeding.cv, verbose=false);

julia> input_fit = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=GenomicBreeding.fit, verbose=false);

julia> input_predict = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, fname_allele_effects_jld2s=["dummy.jld2"], analysis=GenomicBreeding.predict, verbose=false); writejld2(Fit(n=1, l=1), fname="dummy.jld2");

julia> input_gwas = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, analysis=GenomicBreeding.gwas, verbose=false);

julia> inputs_cv = prepareinputs(input_cv); # expect 30 GBInputs = 2 models x 3 traits x (3 populations + 1 bulk + 1 across pops)

julia> inputs_fit = prepareinputs(input_fit); # expect 24 GBInputs = 2 models x 3 traits x (3 populations + 1 bulk)

julia> inputs_predict = prepareinputs(input_predict); # expect 1 GBInput = 1 dummy Fit struct

julia> inputs_gwas = prepareinputs(input_gwas); # expect 24 GBInputs = 1 models x 3 traits x (3 populations + 1 bulk)

julia> length(inputs_cv) == 30
true

julia> length(inputs_fit) == 24
true

julia> length(inputs_predict) == 1
true

julia> length(inputs_gwas) == 24
true

julia> rm.([fname_geno, fname_pheno, "dummy.jld2"]);
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
    _genomes, phenomes, traits_to_skip, populations_to_skip = loadgenomesphenomes(input)
    # Define the model/s to use depending on the type of analysis requested
    models, traits = if input.analysis ∈ [cv, fit]
        input.models, phenomes.traits
    elseif input.analysis ∈ [predict]
        input.fname_allele_effects_jld2s, [nothing]
    elseif input.analysis ∈ [gwas]
        input.gwas_models, phenomes.traits
    else
        # Should be alreat checked in loadgenomesphenomes(...) which calls checkinputs(..)
        valid_analysis_functions = [cv, fit, predict, gwas]
        throw(
            ArgumentError(
                "Analysis: `" *
                string(input.analysis) *
                "` invalid. Please choose from:\n\t‣ " *
                join(string.(valid_analysis_functions), "\n\t‣ "),
            ),
        )
    end
    # Count the number of models, traits, and populations
    m = length(models)
    t = length(traits)
    p = length(unique(phenomes.populations))
    # Extract the populations
    populations = if (p > 1) && (input.analysis ∈ [cv])
        vcat("BULK_CV", "ACROSS_POP_CV", sort(unique(phenomes.populations)))
    elseif (p > 1) && (input.analysis ∈ [fit, gwas])
        vcat("BULK_CV", sort(unique(phenomes.populations)))
    elseif input.analysis ∈ [predict]
        [nothing]
    else
        sort(unique(phenomes.populations))
    end
    # Prepare the GBInputs
    inputs::Vector{GBInput} = []
    # Define the GBInputs for Slurm job arrays or straightforward single computer jobs
    for model in models
        # model = models[1]
        for trait in traits
            # trait = traits[16]
            for population in populations
                # population = populations[4]
                if sum((traits_to_skip .== trait) .&& (populations_to_skip .== population)) > 0
                    # Skip trait-population combinations if we do not have enough entries to perform cross-validation or fit models with confidence
                    continue
                end
                bulk_cv, pops_i = if (population == "BULK_CV") || isnothing(population)
                    true, nothing
                elseif population == "ACROSS_POP_CV"
                    false, populations[3:end]
                else
                    false, [population]
                end
                input_i = clone(input)
                input_i.bulk_cv = bulk_cv
                input_i.populations = pops_i
                input_i.traits = isnothing(trait) ? nothing : [trait]
                if input.analysis ∈ [predict]
                    # Define the model as the filename of the tab-delimited allele effects table
                    # Also set the phenotype file as empty so that we don't merge with a phenomes struct with incomplete correspondence with the genomes struct
                    # Additionally, set the minimum allele frequency to zero so that we don't omit any loci which may have effects in th Fit structs
                    input_i.fname_allele_effects_jld2s = [model]
                    input_i.fname_pheno = ""
                    input_i.maf = 0.0
                else
                    input_i.models = [model]
                end
                push!(inputs, input_i)
            end
        end
    end
    # Output
    inputs
end

"""
    prepareoutprefixandoutdir(input::GBInput)::String

Prepare the output directory and sanitize the output filename prefix for genomic breeding analysis results.

This function performs two main tasks:
1. Creates the output directory if it doesn't exist
2. Sanitizes the output filename prefix by replacing problematic characters with underscores

# Arguments
- `input::GBInput`: Input configuration containing the output file prefix and analysis type

# Returns
- `String`: Sanitized output file prefix path

# Details
Problematic characters that are replaced include spaces, newlines, tabs, parentheses, 
and special characters (`&|:=+*%@!`). The function also ensures the prefix ends with 
the analysis type and a hyphen.

# Throws
- `ArgumentError`: If unable to create the output directory

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase)
julia> input = GBInput(fname_geno="some_dir/fname_geno.jld2", fname_pheno="some_dir/fname_pheno.jld2", fname_out_prefix="GBOutput/some@!_%&prefix", verbose=false);

julia> fname_out_prefix = prepareoutprefixandoutdir(input)
"GBOutput/some_____prefix-cv-"

julia> rm(dirname(fname_out_prefix), recursive=true);
```
"""
function prepareoutprefixandoutdir(input::GBInput)::String
    # Instantiate the output directory if it does not exist
    if !isdir(dirname(input.fname_out_prefix))
        try
            mkdir(dirname(input.fname_out_prefix))
        catch
            throw(ArgumentError("Error creating the output directory: `" * dirname(input.fname_out_prefix) * "`"))
        end
    end
    # Make sure we do not have any problematic symbols in the output filenames
    directory_name = dirname(input.fname_out_prefix)
    prefix_name = if input.fname_out_prefix[end] == '-'
        string(basename(input.fname_out_prefix), input.analysis, "-")
    else
        string(basename(input.fname_out_prefix), "-", input.analysis, "-")
    end
    problematic_strings::Vector{String} = [" ", "\n", "\t", "(", ")", "&", "|", ":", "=", "+", "*", "%", "@", "!"]
    for s in problematic_strings
        prefix_name = replace(prefix_name, s => "_")
    end
    fname_out_prefix = if input.fname_out_prefix != joinpath(directory_name, prefix_name)
        joinpath(directory_name, prefix_name)
    else
        input.fname_out_prefix
    end
    fname_out_prefix
end
