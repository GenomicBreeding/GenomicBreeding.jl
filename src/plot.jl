"""
    plot(;
        input::GBInput,
        skip_genomes::Bool = false,
        skip_phenomes::Bool = false,
        skip_cvs::Bool = false,
        format::String = "svg",
        plot_size::Tuple{Int64,Int64} = (600, 450),
        overwrite::Bool = false
    )::String

Generate and save visualization plots for genomic, phenomic, and cross-validation data.

# Arguments
- `input::GBInput`: Input configuration containing file paths and settings
- `skip_genomes::Bool`: If true, skip generating genome-related plots
- `skip_phenomes::Bool`: If true, skip generating phenome-related plots
- `skip_cvs::Bool`: If true, skip generating cross-validation plots
- `format::String`: Output file format for plots (e.g., "svg", "png")
- `plot_size::Tuple{Int64,Int64}`: Dimensions of output plots in pixels (width, height)
- `overwrite::Bool`: If true, overwrite existing plot files

# Returns
- `String`: Path to the output directory containing generated plots

# Plot Types
## Genomes and Phenomes
- Distribution plots
- Violin plots
- Correlation heatmaps
- Tree plots
- PCA biplots

## Cross-validation
- Bar plots
- Box plots

# Output Structure
Creates a directory structure under `input.fname_out_prefix/plots/` with subdirectories:
- `genomes/`: Genome-related visualizations
- `phenomes/`: Phenome-related visualizations
- `cvs/`: Cross-validation visualizations

# Example
```jldoctest; setup = :(using GenomicBreedingCore, GenomicBreedingIO, GenomicBreeding, StatsBase, DataFrames)
julia> genomes = GenomicBreedingCore.simulategenomes(n=300, l=1_000, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);

julia> trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);

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
    skip_cvs::Bool = false,
    format::String = "svg",
    plot_size::Tuple{Int64,Int64} = (600, 450),
    overwrite::Bool = false,
)::String
    # # genomes = GenomicBreedingCore.simulategenomes(n=300, verbose=false); genomes.populations = StatsBase.sample(string.("pop_", 1:3), length(genomes.entries), replace=true);
    # # trials, _ = GenomicBreedingCore.simulatetrials(genomes=genomes, n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=1, verbose=false);
    # # phenomes = extractphenomes(trials)
    # # fname_geno = try writedelimited(genomes, fname="test-geno.tsv"); catch; rm("test-geno.tsv"); writedelimited(genomes, fname="test-geno.tsv"); end;
    # # fname_pheno = try writedelimited(phenomes, fname="test-pheno.tsv"); catch; rm("test-pheno.tsv"); writedelimited(phenomes, fname="test-pheno.tsv"); end;
    # # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # fname_geno = "test-geno.tsv"; fname_pheno = "test-pheno.tsv"; outdir = "GBOutput"
    # input=GBInput(fname_geno=fname_geno, fname_pheno=fname_pheno, SLURM_cpus_per_task=6, SLURM_mem_G=5)
    # input.fname_out_prefix = replace(input.fname_out_prefix, dirname(input.fname_out_prefix) => outdir)
    # skip_genomes = false; skip_phenomes = false; skip_cvs = false;
    # format = "svg"
    # plot_size = (600, 450); overwrite = true;
    # Load genomes, phenomes, and cvs
    genomes, phenomes = if !skip_genomes || !skip_phenomes
        loadgenomesphenomes(input)
    else
        nothing, nothing
    end
    cvs = if !skip_cvs
        try
            loadcvs(input)
        catch
            nothing
        end
    else
        nothing
    end
    # Define output directory of the plots
    directory_name = if dirname(input.fname_out_prefix) == ""
        pwd()
    else
        dirname(input.fname_out_prefix)
    end
    if !isdir(directory_name)
        try
            mkdir(directory_name)
        catch
            throw(ArgumentError("Error creating the output directory: `" * directory_name * "`"))
        end
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
    if skip_cvs
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
            try
                plots = GenomicBreedingPlots.plot(plot_type, genomes, plot_size = plot_size)
                append!(
                    fnames,
                    saveplots(
                        plots,
                        format = format,
                        prefix = joinpath(plot_outdir, "genomes", string(plot_type)),
                        overwrite = overwrite,
                    ),
                )
            catch
                continue
            end
        end
    end
    # Phenomes
    if !skip_phenomes
        for plot_type in plot_types
            if input.verbose
                println(string("Phenomes: ", plot_type))
            end
            try
                plots = GenomicBreedingPlots.plot(plot_type, phenomes, plot_size = plot_size)
                append!(
                    fnames,
                    saveplots(
                        plots,
                        format = format,
                        prefix = joinpath(plot_outdir, "phenomes", string(plot_type)),
                        overwrite = overwrite,
                    ),
                )
            catch
                continue
            end
        end
    end
    # CVs
    if !skip_cvs && !isnothing(cvs)
        if overwrite
            try
                rm(joinpath(plot_outdir, "cvs"), force = true, recursive = true)
                mkdir(joinpath(plot_outdir, "cvs"))
            catch
                throw(ArgumentError("Error overwriting the output directory: `" * joinpath(plot_outdir, "cvs") * "`"))
            end
        end
        for plot_type in [BarPlots, BoxPlots]
            if input.verbose
                println(string("Vector{CV}: ", plot_type))
            end
            try
                plots = GenomicBreedingPlots.plot(plot_type, cvs, plot_size = plot_size)
                append!(
                    fnames,
                    saveplots(
                        plots,
                        format = format,
                        prefix = joinpath(plot_outdir, "cvs", string(plot_type)),
                        overwrite = overwrite,
                    ),
                )
            catch
                continue
            end
        end
    end
    # Output directory
    if input.verbose
        if length(fnames) > 0
            println(string("Please find output plots in: `", plot_outdir, "`"))
            idx_main_plots = findall(
                .!isnothing.(match.(Regex("PCBiPlots"), fnames)) .||
                .!isnothing.(
                    match.(Regex("BarPlots-Within_population.x.trait.y.cor.z.validation_population.subset."), fnames)
                ),
            )
            println("Please find the following main plots:\n\t‣ " * join(fnames[idx_main_plots], "\n\t‣ "))
        else
            println(string("No plots have been generated. The Slurm jobs may still be running."))
        end
    end
    plot_outdir
end
