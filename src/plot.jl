"""
    plot(;
        input::GBInput,
        skip_genomes::Bool = false,
        skip_phenomes::Bool = false,
        format::String = "svg",
        plot_size::Tuple{Int64,Int64} = (600, 450),
        overwrite::Bool = false,
    )::String

Plot genomes, phenomes, and CVs, if present.

# Example
```jldoctest; setup = :(using GBCore, GBIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GBCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

julia> phenomes = extractphenomes(trials);

julia> fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;

julia> fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    
julia> input = GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5, fname_out_prefix="GBOutput/test-", verbose=false);

julia> GenomicBreeding.plot(input=input, format="png", plot_size = (700, 525))
"GBOutput/plots"

julia> GenomicBreeding.plot(input=input, format="png", plot_size = (700, 525), overwrite=true, skip_genomes=true)
"GBOutput/plots"
```
"""
function plot(;
    input::GBInput,
    skip_genomes::Bool = false,
    skip_phenomes::Bool = false,
    format::String = "svg",
    plot_size::Tuple{Int64,Int64} = (600, 450),
    overwrite::Bool = false,
)::String
    # # genomes = GBCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # # trials, _ = GBCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # # phenomes = extractphenomes(trials)
    # # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # fname_geno = "test-geno.tsv"; fname_pheno = "test-pheno.tsv"; outdir = "GBOutput"
    # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # input.fname_out_prefix = replace(input.fname_out_prefix, dirname(input.fname_out_prefix) => outdir)
    # format = "svg"
    # Load genomes and phenomes
    genomes, phenomes = if !skip_genomes || !skip_phenomes
        load(input)
    else
        nothing, nothing
    end
    # Define output directory of the CVs
    directory_name = if dirname(input.fname_out_prefix) == ""
        pwd()
    else
        dirname(input.fname_out_prefix)
    end
    fnames_cvs = if isdir(directory_name)
        files = readdir(directory_name)
        idx = findall(.!isnothing.(match.(Regex("jld2\$"), files)))
        fnames_cvs = if length(idx) > 0
            joinpath.(directory_name, files[idx])
        else
            nothing
        end
        fnames_cvs
    else
        try
            mkdir(directory_name)
        catch
            throw(ArgumentError("Error creating the output directory: `" * directory_name * "`"))
        end
        nothing
    end
    # Prepare the main directory for the output plots
    plot_outdir = joinpath(directory_name, "plots")
    if !isdir(plot_outdir)
        try
            mkdir(plot_outdir)
        catch
            throw(ArgumentError("Error creating the output directory: `" * plot_outdir * "`"))
        end
    end
    subdir_names = ["genomes", "phenomes", "cvs"]
    if skip_genomes
        subdir_names = subdir_names[subdir_names.!="genomes"]
    end
    if skip_phenomes
        subdir_names = subdir_names[subdir_names.!="phenomes"]
    end
    if isnothing(fnames_cvs)
        subdir_names = subdir_names[subdir_names.!="cvs"]
    end
    for subdir_name in subdir_names
        # subdir_name = subdir_names[1]
        plot_outdir_subdir = joinpath(plot_outdir, subdir_name)
        if !isdir(plot_outdir_subdir)
            try
                mkdir(plot_outdir_subdir)
            catch
                throw(ArgumentError("Error creating the output directory: `" * plot_outdir_subdir * "`"))
            end
        end
    end
    # Output plot filenames
    fnames::Vector{String} = []
    # Plot types for Genomes and Phenomes
    plot_types = [DistributionPlots, ViolinPlots, CorHeatPlots, TreePlots, PCBiPlots]
    # Genomes
    if !skip_genomes
        for plot_type in plot_types
            # plot_type = PCBiPlots
            if input.verbose
                println(string("Genomes: ", plot_type))
            end
            plots = GBPlots.plot(plot_type, genomes, plot_size = plot_size)
            append!(
                fnames,
                saveplots(
                    plots,
                    format = format,
                    prefix = joinpath(plot_outdir, "genomes", string(plot_type)),
                    overwrite = overwrite,
                ),
            )
        end
    end
    # Phenomes
    if !skip_phenomes
        for plot_type in plot_types
            if input.verbose
                println(string("Phenomes: ", plot_type))
            end
            plots = GBPlots.plot(plot_type, phenomes, plot_size = plot_size)
            append!(
                fnames,
                saveplots(
                    plots,
                    format = format,
                    prefix = joinpath(plot_outdir, "phenomes", string(plot_type)),
                    overwrite = overwrite,
                ),
            )
        end
    end
    # CVs
    if !isnothing(fnames_cvs)
        if overwrite
            try
                rm(joinpath(plot_outdir, "cvs"), force = true, recursive = true)
                mkdir(joinpath(plot_outdir, "cvs"))
            catch
                throw(ArgumentError("Error overwriting the output directory: `" * joinpath(plot_outdir, "cvs") * "`"))
            end
        end
        cvs = Vector{CV}(undef, length(fnames_cvs))
        for (i, fname) in enumerate(fnames_cvs)
            # i = 1; fname = fnames_cvs[i];
            cvs[i] = readjld2(CV, fname = fname)
        end
        for plot_type in [BarPlots, BoxPlots]
            if input.verbose
                println(string("Vector{CV}: ", plot_type))
            end
            plots = GBPlots.plot(plot_type, cvs, plot_size = plot_size)
            append!(
                fnames,
                saveplots(
                    plots,
                    format = format,
                    prefix = joinpath(plot_outdir, "cvs", string(plot_type)),
                    overwrite = overwrite,
                ),
            )
        end
    end
    # Output directory
    if input.verbose
        if length(fnames) > 0
            println(string("Please find output plots in: `", plot_outdir, "`"))
            println(fnames)
        else
            println(string("No plots have been generated. The Slurm jobs may still be running."))
        end
    end
    plot_outdir
end
