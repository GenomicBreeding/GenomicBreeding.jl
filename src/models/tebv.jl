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
function trialsmodelsfomulae!(df::DataFrame; trait::String, max_levels::Int64 = 100)::Tuple{Vector{String}, Vector{Int64}}
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
    vars::Array{Array{String,1},1} = [[x] for x in explanatory_variables]
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
                if (spatial_error_column != "") && (countlevels(df; column_names = [interaction_column, spatial_error_column]) <= max_levels)
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
                    push!(n_levels, countlevels(df; column_names = vcat(["entries"], fixed, interaction_column, spatial_error_column)))
                end
            else
                if countlevels(df; column_names = [interaction_column]) <= max_levels
                    push!(formulae, string(trait, " ~ 1 + (0 + ", interaction_column, "|entries)"))
                    push!(n_levels, countlevels(df; column_names = ["entries", interaction_column]))
                end
                if (spatial_error_column != "") && (countlevels(df; column_names = [interaction_column, spatial_error_column]) <= max_levels)
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
                    push!(n_levels, countlevels(df; column_names = ["entries", interaction_column, spatial_error_column]))
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
function analyse(trials::Trials; max_levels::Int64 = 100, max_time_per_model::Int64 = 60)::Bool
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes()); max_levels::Int64=100; max_time_per_model::Int64=60;
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes(n=5), n_years=2, n_seasons=2, n_harvests=1, n_sites=2, n_replications=10); max_levels::Int64=100; max_time_per_model::Int64=60;
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes(n=5), n_years=1, n_seasons=1, n_harvests=1, n_sites=1, n_replications=10); max_levels::Int64=100; max_time_per_model::Int64=60;
    # Tabularise
    df::DataFrame = tabularise(trials)
    # Number of entries whose BLUEs or BLUPs we wish to extract
    entries::Vector{String} = sort(unique(df[!, "entries"]))
    n::Int64 = length(entries)
    # Iterate across each trait
    for trait in unique(trials.traits)
        # trait = unique(trials.traits)[1]
        # Define the formulae
        formulae, n_levels = trialsmodelsfomulae!(df; trait = trait, max_levels = max_levels)
        # Fit the models
        models::Vector{LinearMixedModel{Float64}} = Vector{LinearMixedModel{Float64}}(undef, length(formulae))
        deviances::Vector{Float64} = fill(Inf64, length(formulae))
        # Parallel fitting of models with total number of levels at or below `(1.5 * max_levels)`
        @time Threads.@threads for i = findall(n_levels .<= (1.5 * max_levels))
            f = @eval(@string2formula $(formulae[i]))
            model = MixedModel(f, df)
            model.optsum.REML = true
            model.optsum.maxtime = max_time_per_model
            @time fit!(model)
            models[i] = model
            deviances[i] = deviance(model)
        end
        # Iterative fitting of models with total number of levels above `(1.5 * max_levels)` to minimise OOM errors from occuring
        @time for i = findall((n_levels .<= (1.5 * max_levels)) .!= true)
            f = @eval(@string2formula $(formulae[i]))
            model = MixedModel(f, df)
            model.optsum.REML = true
            model.optsum.maxtime = max_time_per_model
            @time fit!(model)
            models[i] = model
            deviances[i] = deviance(model)
            println(string("[", i, "]: ", f))
        end
        model = models[argmin(deviances)]
        # model = models[sample(1:length(deviances), 1)[1]]
        deviances[argmin(deviances)]
        formulae[argmin(deviances)]
        # Model stats
        plt = UnicodePlots.lineplot(deviances)
        UnicodePlots.hline!(plt, deviances[argmin(deviances)]; name = "minimum")
        model.optsum.returnvalue
        # "y = Xβ + Zb"
        # X = model.X
        # β = model.β
        # Z = only(model.reterms)
        # b = model.b
        # Σ = varest(model) .* (only(model.λ) * only(model.λ)')
        # σ = model.σ


        cors = []
        for i in eachindex(models)
            println(string("[", i, "]: ", formulae[i]))
            model = models[i]
            # Fixed effects
            df_BLUEs = DataFrame(coeftable(model))
            idx_row = findall(match.(r"entries: ", df_BLUEs[!, "Name"]) .!= nothing)
            if length(idx_row) == (n - 1)
                df_BLUEs[idx_row, "Name"] = replace.(df_BLUEs[idx_row, "Name"], "entries: " => "")
                df_BLUEs[1, "Name"] = setdiff(entries, df_BLUEs[idx_row, "Name"])[1]
                df_BLUEs[2:end, "Coef."] .+= df_BLUEs[1, "Coef."]
                DataFrame(; entries = entries)
            end
            # BLUPs
            df_BLUPs = DataFrame(only(raneftables(model)))
            # Find the columns referring to the nested entry effects, and excluding the spatial residual effects
            idx_col = findall(match.(r"blocks|rows|cols", names(df_BLUPs)) .== nothing)
            df_BLUPs = df_BLUPs[:, idx_col]
            g_hat = sum.(eachrow(df_BLUPs[:, 2:end]))

            if length(idx_row) == (n - 1)
                g_hat = df_BLUEs[!, "Coef."]
            end

            UnicodePlots.histogram(g_hat)

            additive_genetic = simulated_effects[1].additive_genetic
            for i = 2:length(simulated_effects)
                additive_genetic += simulated_effects[i].additive_genetic
            end
            additive_genetic ./= length(simulated_effects)
            UnicodePlots.histogram(additive_genetic)
            
            y = g_hat
            X = hcat(ones(length(additive_genetic)), additive_genetic)
            b_hat = X \ y
            p = 1_000
            X_new = hcat(ones(p), collect(range(minimum(additive_genetic), maximum(additive_genetic), length=p)))
            y_hat = X_new * b_hat

            plt = UnicodePlots.scatterplot(additive_genetic, g_hat)
            UnicodePlots.lineplot!(plt, X_new[:,2], y_hat)
            append!(cors, cor(g_hat, additive_genetic))
        end

        UnicodePlots.scatterplot(deviances, cors)


    end

    return false
end
