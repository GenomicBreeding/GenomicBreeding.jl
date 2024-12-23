"""
    Base.hash(x::TEBV, h::UInt)::UInt

Hash a TEBV struct using the traits, formualae and phenomes.
We deliberately excluded the models, df_BLUEs, and df_BLUPs for efficiency.

## Examples
```jldoctest; setup = :(using GenomicBreeding, MixedModels, DataFrames)
julia> tebv = TEBV(traits=[""], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1,t=1)]);

julia> typeof(hash(tebv))
UInt64
```
"""
function Base.hash(x::TEBV, h::UInt)::UInt
    # hash(TEBV, hash(x.traits, hash(x.formulae, hash(x.models, hash(x.df_BLUEs, hash(x.df_BLUPs, hash(x.phenomes, h)))))))
    hash(TEBV, hash(x.traits, hash(x.formulae, hash(x.phenomes, h))))
end


"""
    Base.:(==)(x::TEBV, y::TEBV)::Bool

Equality of TEBV structs using the hash function defined for TEBV structs.

## Examples
```jldoctest; setup = :(using GenomicBreeding, MixedModels, DataFrames)
julia> tebv_1 = TEBV(traits=[""], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1,t=1)]);

julia> tebv_2 = TEBV(traits=[""], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1,t=1)]);

julia> tebv_3 = TEBV(traits=["SOMETHING_ELSE"], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1,t=1)]);

julia> tebv_1 == tebv_2
true

julia> tebv_1 == tebv_3
false
```
"""
function Base.:(==)(x::TEBV, y::TEBV)::Bool
    hash(x) == hash(y)
end


"""
    checkdims(y::TEBV)::Bool

Check dimension compatibility of the fields of the TEBV struct

## Examples
```jldoctest; setup = :(using GenomicBreeding, MixedModels, DataFrames)
julia> tebv = TEBV(traits=[""], formulae=[""], models=[MixedModel(@formula(y~1+(1|x)), DataFrame(y=1, x=1))], df_BLUEs=[DataFrame(x=1)], df_BLUPs=[DataFrame(x=1)], phenomes=[Phenomes(n=1,t=1)]);

julia> checkdims(tebv)
true
```
"""
function checkdims(tebv::TEBV)::Bool
    t = length(tebv.traits)
    if (t != length(tebv.formulae)) ||
       (t != length(unique(tebv.models))) ||
       (t != length(unique(tebv.df_BLUEs))) ||
       (t != length(unique(tebv.df_BLUPs))) ||
       (t != length(unique(tebv.phenomes)))
        return false
    end
    if !isa(tebv.traits, Vector{String}) ||
       !isa(tebv.formulae, Vector{String}) ||
       !isa(tebv.models, Vector{LinearMixedModel{Float64}}) ||
       !isa(tebv.df_BLUEs, Vector{DataFrame}) ||
       !isa(tebv.df_BLUPs, Vector{DataFrame}) ||
       !isa(tebv.phenomes, Vector{Phenomes})
        return false
    end
    return true
end


"""
    countlevels(df::DataFrame; column_names::Vector{String})::Int64

Count the total number of factor levels across the column names provided.
"""
function countlevels(df::DataFrame; column_names::Vector{String})::Int64
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes()); df::DataFrame = tabularise(trials); column_names::Vector{String} = ["years", "entries"]
    if length(names(df) ∩ column_names) != length(column_names)
        throw(ArgumentError("The supplied column names are not found in the dataframe."))
    end
    m::Int64 = 0
    for column_name in column_names
        m += length(unique(df[!, column_name]))
    end
    return m
end


"""
    string2formula(x)

Macro to `Meta.parse` a string into a formula.
"""
macro string2formula(x)
    @eval(@formula($(Meta.parse(x))))
end


"""
    trialsmodelsfomulae!(df::DataFrame; trait::String, max_levels::Int64 = 100)::Vector{String}

Define formulae for the mixed models to fit on the tabularised `Trials` struct.
    - appends interaction effects intto `df`
    - returns:
        + a vector of formulae as strings
        + a vector of the total number of non-entry factor levels

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials, _simulated_effects = simulatetrials(genomes = simulategenomes(verbose=false), verbose=false);

julia> df = tabularise(trials);

julia> size(df)
(12800, 14)

julia> formulae, n_levels = trialsmodelsfomulae!(df, trait="trait_1");

julia> size(df)
(12800, 134)

julia> length(formulae)
76

julia> sum(n_levels .== sort(n_levels))
76
```
"""
function trialsmodelsfomulae!(
    df::DataFrame;
    trait::String,
    max_levels::Int64 = 100,
)::Tuple{Vector{String},Vector{Int64}}
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes()); df::DataFrame = tabularise(trials); trait::String=trials.traits[1]; max_levels::Int64=100;
    # Define the totality of all expected variables
    nester_variables = ["years", "seasons", "harvests", "sites"]
    spatial_variables = ["blocks", "rows", "cols"]
    target_variables = ["entries", "populations"]
    residual_variable = ["replications"]
    # Extract the available (non-fixed) subset of all the variables listed above
    explanatory_variables = filter(
        x -> length(unique(df[!, x])) > 1,
        vcat(nester_variables, spatial_variables, target_variables, residual_variable),
    )
    # Warn if unreplicated
    if sum(explanatory_variables .== "replications") == 0
        @warn "Unreplicated trial data!"
    end
    # We expect to extract BLUEs/BLUPs of the entries
    if sum(explanatory_variables .== "entries") == 0
        throw(ErrorException("Only one entry in the entire trial!"))
    end
    # Exclude replications (residuals estimator), entries (what we're most interested in), and population (entry-specific and can be used instead of entries) in the formula
    explanatory_variables = filter(x -> x != "replications", explanatory_variables)
    explanatory_variables = filter(x -> x != "entries", explanatory_variables)
    explanatory_variables = filter(x -> x != "populations", explanatory_variables)
    # Permute
    vars::Array{Vector{String},1} = [[x] for x in explanatory_variables]
    for nv in explanatory_variables
        tmp_v = []
        for v in vars
            v = copy(v)
            push!(v, nv)
            sort!(v)
            unique!(v)
            push!(tmp_v, v)
        end
        sort!(tmp_v)
        unique!(tmp_v)
        append!(vars, tmp_v)
    end
    sort!(vars)
    unique!(vars)
    # Define models where we nest entries, blocks, rows and cols
    formulae::Vector{String} = []
    n_levels::Vector{Int64} = []
    for var in vars
        # var = vars[7]
        # Sort by decreasing nester
        spatials = var[[findall(var .== v)[1] for v in spatial_variables ∩ var]]
        nesters = var[[findall(var .== v)[1] for v in nester_variables ∩ var]]
        # Define the spatial errors nested within the nester variables
        spatial_error_column::String = ""
        if length(spatials) > 0
            spatial_error_column = join(vcat(nesters, spatials), "_x_")
            try
                df[!, spatial_error_column]
            catch
                nesters_spatials = vcat(nesters[2:end], spatials)
                df[!, spatial_error_column] = df[!, nesters_spatials[1]]
                for v in nesters_spatials[2:end]
                    df[!, spatial_error_column] = string.(df[!, spatial_error_column], "_x_", df[!, v])
                end
            end
        end

        # Single environment model
        if length(nesters) == 0
            if spatial_error_column != ""
                if countlevels(df; column_names = [spatial_error_column]) <= max_levels
                    # Fixed entries
                    push!(formulae, string(trait, " ~ 1 + entries + (0 + ", spatial_error_column, "|entries)"))
                    push!(n_levels, countlevels(df; column_names = ["entries", spatial_error_column]))
                    # Random entries
                    push!(formulae, string(trait, " ~ 1 + (1|entries) + (0 + ", spatial_error_column, "|entries)"))
                    push!(n_levels, countlevels(df; column_names = ["entries", spatial_error_column]))
                end
            else
                # Fixed entry effects
                # ... BLUEs would simply be means.
                # Random entries
                push!(formulae, string(trait, " ~ 1 + (0 + 1|entries)"))
                push!(n_levels, countlevels(df; column_names = ["entries"]))
            end
        end

        # Multi-environemnt models
        for i in eachindex(nesters)
            # i = 2
            # Divide the variables into fixed effects and random interaction effects
            fixed::Vector{String} = nesters[1:(end-i)]
            random_interaction::Vector{String} = nesters[(length(nesters)-(i-1)):end]
            if length(random_interaction) > 1
                interaction_column::String = join(random_interaction, "_x_")
                try
                    df[!, interaction_column]
                catch
                    df[!, interaction_column] = df[!, random_interaction[1]]
                    for v in random_interaction[2:end]
                        df[!, interaction_column] = string.(df[!, interaction_column], "_x_", df[!, v])
                    end
                end
                if length(unique(df[!, interaction_column])) > max_levels
                    continue
                end
            else
                interaction_column = random_interaction[1]
            end
            # Append the model
            if length(fixed) > 0
                if countlevels(df; column_names = [interaction_column]) <= max_levels
                    push!(
                        formulae,
                        string(trait, " ~ 1 + ", join(fixed, " + "), " + (0 + ", interaction_column, "|entries)"),
                    )
                    push!(n_levels, countlevels(df; column_names = vcat(["entries"], fixed, interaction_column)))
                end
                if (spatial_error_column != "") &&
                   (countlevels(df; column_names = [interaction_column, spatial_error_column]) <= max_levels)
                    push!(
                        formulae,
                        string(
                            trait,
                            " ~ 1 + ",
                            join(fixed, " + "),
                            " + (0 + ",
                            interaction_column,
                            "|entries)",
                            " + (0 + ",
                            spatial_error_column,
                            "|entries)",
                        ),
                    )
                    push!(
                        n_levels,
                        countlevels(
                            df;
                            column_names = vcat(["entries"], fixed, interaction_column, spatial_error_column),
                        ),
                    )
                end
            else
                if countlevels(df; column_names = [interaction_column]) <= max_levels
                    push!(formulae, string(trait, " ~ 1 + (0 + ", interaction_column, "|entries)"))
                    push!(n_levels, countlevels(df; column_names = ["entries", interaction_column]))
                end
                if (spatial_error_column != "") &&
                   (countlevels(df; column_names = [interaction_column, spatial_error_column]) <= max_levels)
                    push!(
                        formulae,
                        string(
                            trait,
                            " ~ 1 + (0 + ",
                            interaction_column,
                            "|entries)",
                            " + (0 + ",
                            spatial_error_column,
                            "|entries)",
                        ),
                    )
                    push!(
                        n_levels,
                        countlevels(df; column_names = ["entries", interaction_column, spatial_error_column]),
                    )
                end
            end
        end
    end
    # Keep only unique formulae
    idx = unique(i -> formulae[i], eachindex(formulae))
    formulae = formulae[idx]
    n_levels = n_levels[idx]
    # Sort formulae by complexity-ish
    idx = sortperm(n_levels)
    formulae = formulae[idx]
    n_levels = n_levels[idx]
    # Remove models with blocks if rows and columns exist
    if length(unique(vcat(vars...)) ∩ ["rows", "cols"]) == 2
        idx = findall(match.(r"blocks", formulae) .== nothing)
        formulae = formulae[idx]
        n_levels = n_levels[idx]
    end
    # Remove models with redundant nesters in the main and residual terms
    idx = findall([
        (sum(match.(r"years_x_seasons_x_harvests_x_sites", x) .!= nothing) < 2) &&
        (sum(match.(r"years_x_seasons_x_sites", x) .!= nothing) < 2) for x in split.(formulae, " + (")
    ])
    formulae = formulae[idx]
    n_levels = n_levels[idx]
    # Output in addition to the mutated `df`
    (formulae, n_levels)
end

# NOTES: 
# (1) Avoid over-parameterisation we'll have enough of that with the genomic prediction models
# (2) We will fit mixed models with unstructure variance-covariance matrix of the random effects
# (3) We prefer REML over ML
# (4) We will compare BLUEs vs BLUPs of entries
"""
# Analyse trials

## Arguments
- `trials`: Trials struct 
- `max_levels`: maximum number of non-entry factor levels to include in the linear mixed models (default = 100)
- `max_time_per_model`: maximum time in seconds for fitting each linear mixed model (default = 60)
- `verbose`: Show trials analysis progress bar? (default = true)

## Outputs
- `TEBV` struct containing the trait names, the best fitting formulae, models, BLUEs, and BLUPs for each trait

## Examples
```jldoctest; setup = :(using GenomicBreeding; using StatsBase)
julia> trials, _simulated_effects = simulatetrials(genomes = simulategenomes(n=10, verbose=false), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10, verbose=false);

julia> tebv = analyse(trials, max_levels=50, verbose=false);

julia> tebv.traits
3-element Vector{String}:
 "trait_1"
 "trait_2"
 "trait_3"

julia> checkdims(tebv)
true
```
"""
function analyse(trials::Trials; max_levels::Int64 = 100, max_time_per_model::Int64 = 60, verbose::Bool = true)::TEBV
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes()); max_levels::Int64=100; max_time_per_model::Int64=60; verbose::Bool = true;
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes(n=5), n_years=2, n_seasons=2, n_harvests=1, n_sites=2, n_replications=10); max_levels::Int64=100; max_time_per_model::Int64=60; verbose::Bool = true;
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes(n=5), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10); max_levels::Int64=100; max_time_per_model::Int64=60; verbose::Bool = true;
    # Tabularise
    df::DataFrame = tabularise(trials)
    # Number of entries whose BLUEs or BLUPs we wish to extract
    entries::Vector{String} = sort(unique(df[!, "entries"]))
    n::Int64 = length(entries)
    if verbose
        println("Analysing trial data with:")
        println(string("\t ‣ total observations = ", nrow(df)))
        println(string("\t ‣ entries = ", n))
        println(string("\t ‣ populations = ", length(unique(trials.populations))))
        println(string("\t ‣ traits = ", length(trials.traits)))
        println(string("\t ‣ years = ", length(unique(trials.years))))
        println(string("\t ‣ seasons = ", length(unique(trials.seasons))))
        println(string("\t ‣ harvests = ", length(unique(trials.harvests))))
        println(string("\t ‣ sites = ", length(unique(trials.sites))))
        println(string("\t ‣ replications = ", length(unique(trials.replications))))
        println(string("\t ‣ blocks = ", length(unique(trials.blocks))))
        println(string("\t ‣ rows = ", length(unique(trials.rows))))
        println(string("\t ‣ cols = ", length(unique(trials.cols))))
    end
    # Iterate across each trait
    out_traits::Vector{String} = []
    out_formulae::Vector{String} = []
    out_models::Vector{LinearMixedModel{Float64}} = []
    out_df_BLUEs::Vector{DataFrame} = []
    out_df_BLUPs::Vector{DataFrame} = []
    out_phenomes::Vector{Phenomes} = []
    for trait in unique(trials.traits)
        # trait = unique(trials.traits)[3]
        if verbose
            println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            println(trait)
            println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        end
        # Define the formulae
        formulae, n_levels = trialsmodelsfomulae!(df; trait = trait, max_levels = max_levels)
        idx_parallel_models::Vector{Int64} = findall(n_levels .<= (1.5 * max_levels))
        idx_iterative_models::Vector{Int64} = findall((n_levels .<= (1.5 * max_levels)) .!= true)
        if verbose
            println(string("\t ‣ models to fit = ", length(formulae)))
            println(string("\t ‣ parallel fitting for the ", length(idx_parallel_models), " simplest models"))
            println(string("\t ‣ iterative fitting for the ", length(idx_iterative_models), " most complex models"))
            println(string("\t ‣ simplest model: ", formulae[1]))
            println(string("\t ‣ most complex model: ", formulae[end]))
        end
        # Fit the models
        models::Vector{LinearMixedModel{Float64}} = Vector{LinearMixedModel{Float64}}(undef, length(formulae))
        deviances::Vector{Float64} = fill(Inf64, length(formulae))
        # Parallel fitting of models with total number of levels at or below `(1.5 * max_levels)`
        if verbose
            pb = Progress(length(idx_parallel_models); desc = "Trials analyses | parallel model fitting: ")
        end
        thread_lock::ReentrantLock = ReentrantLock()
        Threads.@threads for i in idx_parallel_models
            # for i in idx_parallel_models
            f = @eval(@string2formula $(formulae[i]))
            model = MixedModel(f, df)
            @lock thread_lock model.optsum.REML = true
            @lock thread_lock model.optsum.maxtime = max_time_per_model
            @lock thread_lock try
                fit!(model, progress = false)
            catch
                try
                    fit!(model, progress = false)
                catch
                    continue
                end
            end
            @lock thread_lock models[i] = model
            @lock thread_lock deviances[i] = deviance(model)
            if verbose
                next!(pb)
            end
        end
        if verbose
            finish!(pb)
        end
        # Iterative fitting of models with total number of levels above `(1.5 * max_levels)` to minimise OOM errors from occuring
        if verbose
            pb = Progress(length(idx_parallel_models); desc = "Trials analyses | iterative model fitting: ")
        end
        for i in idx_iterative_models
            f = @eval(@string2formula $(formulae[i]))
            model = MixedModel(f, df)
            model.optsum.REML = true
            model.optsum.maxtime = max_time_per_model
            try
                fit!(model, progress = false)
            catch
                try
                    fit!(model, progress = false)
                catch
                    continue
                end
            end
            models[i] = model
            deviances[i] = deviance(model)
            if verbose
                next!(pb)
            end
        end
        if verbose
            finish!(pb)
        end
        model = models[argmin(deviances)]
        # model = models[sample(1:length(deviances), 1)[1]]
        # # Model stats
        # plt = UnicodePlots.lineplot(deviances)
        # UnicodePlots.hline!(plt, deviances[argmin(deviances)]; name = "minimum")
        # model.optsum.returnvalue
        # deviances[argmin(deviances)]
        # # "y = Xβ + Zb"
        # # X = model.X
        # # β = model.β
        # # Z = only(model.reterms)
        # # b = model.b
        # # Σ = varest(model) .* (only(model.λ) * only(model.λ)')
        # # σ = model.σ
        # Fixed effects
        df_BLUEs::DataFrame = DataFrame(coeftable(model))
        # Determine the effect of the first entry and add the intercept to the effects of the other entries
        idx_row = findall(match.(r"entries: ", df_BLUEs[!, "Name"]) .!= nothing)
        if length(idx_row) == (n - 1)
            df_BLUEs[idx_row, "Name"] = replace.(df_BLUEs[idx_row, "Name"], "entries: " => "")
            df_BLUEs[1, "Name"] = setdiff(entries, df_BLUEs[idx_row, "Name"])[1]
            df_BLUEs[2:end, "Coef."] .+= df_BLUEs[1, "Coef."]
        end
        # BLUPs
        df_BLUPs::DataFrame = DataFrame(only(raneftables(model)))
        # Find non-spatial interaction effects with the entries
        idx_col::Vector{Int64} = findall(match.(r"blocks|rows|cols", names(df_BLUPs)) .== nothing)
        # Instantiate the Phenomes struct
        t::Int64 = 1
        if length(idx_col) > 1
            t = length(idx_col)
        end
        phenomes::Phenomes = Phenomes(n = length(entries), t = t)
        phenomes.mask .= true
        # If we have BLUPs (multiple columns or "multiple traits")
        if length(idx_col) > 1
            phenomes.entries = df_BLUPs[:, 1]
            phenomes.phenotypes = Matrix(df_BLUPs[:, idx_col[2:end]])
            phenomes.traits = string.(trait, "|", names(df_BLUPs))[idx_col[2:end]]
        end
        # If we have BLUEs (1 column or "1 trait" only)
        if length(idx_row) == (n - 1)
            phenomes.entries = df_BLUEs[!, "Name"]
            phenomes.phenotypes[:, 1] = df_BLUEs[!, "Coef."]
            phenomes.traits[1] = trait
        end
        # Populate populations
        for (i, entry) in enumerate(phenomes.entries)
            idx = findall(trials.entries .== entry)[1]
            phenomes.populations[i] = trials.populations[idx]
        end
        # Update output vectors
        push!(out_traits, trait)
        push!(out_formulae, formulae[argmin(deviances)])
        push!(out_models, model)
        push!(out_df_BLUEs, df_BLUEs)
        push!(out_df_BLUPs, df_BLUPs)
        push!(out_phenomes, phenomes)
    end
    TEBV(
        traits = out_traits,
        formulae = out_formulae,
        models = out_models,
        df_BLUEs = out_df_BLUEs,
        df_BLUPs = out_df_BLUPs,
        phenomes = out_phenomes,
    )
end

function plots(trials::Trials)::Vector{Plots} end