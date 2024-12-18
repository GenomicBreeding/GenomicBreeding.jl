"""
    definetrialsmodelsfomulae!(df::DataFrame; trait::String, max_levels::Int64 = 100)::Vector{String}

Define formulae for the mixed models to fit on the tabularised `Trials` struct.
    - appends interaction effects intto `df`
    - returns a vector of formulae as strings

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials, _simulated_effects = simulatetrials(genomes = simulategenomes(verbose=false), verbose=false);

julia> df = tabularise(trials);

julia> size(df)
(12800, 14)

julia> formulae = definetrialsmodelsfomulae!(df, trait="trait_1");

julia> size(df)
(12800, 134)

julia> length(formulae)
166
```
"""
function definetrialsmodelsfomulae!(
    df::DataFrame;
    trait::String,
    max_levels::Int64 = 100,
)::Vector{String}
    # trials, simulated_effects = simulatetrials(genomes = simulategenomes());
    # df::DataFrame = tabularise(trials); trait::String=trials.traits[1]; max_levels::Int64=100;
    nester_variables = ["years", "seasons", "harvests", "sites"]
    spatial_variables = ["blocks", "rows", "cols"]
    target_variables = ["entries", "populations"]
    residual_variable = ["replications"]
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
    # Find nester variables and permutate
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
    for var in vars
        # var = vars[7]

        # Sort by decreasing nester
        spatials = var[[findall(var .== v)[1] for v in spatial_variables ∩ var]]
        nesters = var[[findall(var .== v)[1] for v in nester_variables ∩ var]]

        # Define the spatial errors nested within the nester variables
        spatial_error_column::String = join(vcat(nesters, spatials), "_x_")
        if spatial_error_column != ""
            try
                df[!, spatial_error_column]
            catch
                nesters_spatials = vcat(nesters[2:end], spatials)
                df[!, spatial_error_column] = df[!, nesters_spatials[1]]
                for v in nesters_spatials[2:end]
                    df[!, spatial_error_column] =
                        string.(df[!, spatial_error_column], "_x_", df[!, v])
                end
            end
        end

        # Single environment model
        if length(nesters) == 0
            if spatial_error_column != ""
                # Fixed entries
                push!(
                    formulae,
                    string(
                        trait,
                        " ~ 1 + entries + (0 + ",
                        spatial_error_column,
                        "|entries)",
                    ),
                )
                # Random entries
                push!(
                    formulae,
                    string(
                        trait,
                        " ~ 1 + (1|entries) + (0 + ",
                        spatial_error_column,
                        "|entries)",
                    ),
                )
            else
                # Fixed entrie effects
                # ... BLUEs would simply be means.
                # Random entries
                push!(formulae, string(trait, " ~ 1 + (0 + 1|entries)"))
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
                        df[!, interaction_column] =
                            string.(df[!, interaction_column], "_x_", df[!, v])
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
                push!(
                    formulae,
                    string(
                        trait,
                        " ~ 1 + ",
                        join(fixed, " + "),
                        " + (0 + ",
                        interaction_column,
                        "|entries)",
                    ),
                )
                if spatial_error_column != ""
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
                end
            else
                # push!(formulae, string(trait, " ~ 1 + ", interaction_column, "|entries"))
                push!(
                    formulae,
                    string(trait, " ~ 1 + (0 + ", interaction_column, "|entries)"),
                )
                if spatial_error_column != ""
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
                end
            end
        end
    end
    sort!(formulae)
    unique!(formulae)
    # Sort formulae by complexity-ish
    formulae = formulae[sortperm(length.(formulae))]
    # Remove models with blocks if rows and columns exist
    if length(unique(vcat(vars...)) ∩ ["rows", "cols"]) == 2
        formulae = formulae[match.(r"blocks", formulae).==nothing]
    end
    # Remove models with redundant nesters in the main and residual terms
    idx = findall(
        [
            (sum(match.(r"years_x_seasons_x_harvests_x_sites", x) .!= nothing) < 2) && 
            (sum(match.(r"years_x_seasons_x_sites", x) .!= nothing) < 2) 
            for x in split.(formulae, " + (")
        ]
    )
    formulae = formulae[idx]
    # Output in addition to the mutated `df`
    formulae
end

macro string2formula(x)
     @eval(@formula($(Meta.parse(x))))
end

# NOTES: 
# (1) Avoid over-parameterisation we'll have enough of that with the genomic prediction models
# (2) We will fit mixed models with unstructure variance-covariance matrix of the random effects
# (3) We prefer REML over ML
# (4) We will compare BLUEs vs BLUPs of entries
function analyse(
    trials::Trials;
    max_levels::Int64 = 100,
    max_time_per_model::Int64 = 60,
)::Bool
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
        formulae = definetrialsmodelsfomulae!(df, trait = trait, max_levels = max_levels)
        # Fit the models
        models::Vector{LinearMixedModel{Float64}} =
            Vector{LinearMixedModel{Float64}}(undef, length(formulae))
        deviances::Vector{Float64} = fill(Inf64, length(formulae))
        if length(formulae) > 15
            # Start with the 10 least saturated models
            @time Threads.@threads for i = 1:minimum([10, length(formulae)])
                f = @eval(@string2formula $(formulae[i]))
                model = MixedModel(f, df)
                model.optsum.REML = true
                model.optsum.maxtime = max_time_per_model
                @time fit!(model)
                models[i] = model
                deviances[i] = deviance(model)
            end
            # Then proceed with the 5 most saturated models, 
            # but instead of fitting in parallel, we do this iteratively to minimise OOM errors
            @time for i = (length(formulae)-4):length(formulae)
                f = @eval(@string2formula $(formulae[i]))
                model = MixedModel(f, df)
                model.optsum.REML = true
                model.optsum.maxtime = max_time_per_model
                @time fit!(model)
                models[i] = model
                deviances[i] = deviance(model)
            end
        else
            # Fit all models if we have at most 15 formulae
            @time Threads.@threads for i in eachindex(formulae)
                # i = 2
                f = @eval(@string2formula $(formulae[i]))
                model = MixedModel(f, df)
                model.optsum.maxtime = max_time_per_model
                @time fit!(model, REML = true)
                models[i] = model
                deviances[i] = deviance(model)
            end
        end
        model = models[argmin(deviances)]
        # Model stats
        plt = UnicodePlots.lineplot(deviances)
        UnicodePlots.hline!(plt, deviances[argmin(deviances)], name = "minimum")
        model.optsum.returnvalue
        # "y = Xβ + Zb"
        # X = model.X
        # β = model.β
        # Z = only(model.reterms)
        # b = model.b
        # Σ = varest(model) .* (only(model.λ) * only(model.λ)')
        # σ = model.σ

        # Fixed effects
        df_BLUEs = DataFrame(coeftable(model))
        idx = findall(match.(r"entries: ", df_BLUEs[!, "Name"]) .!= nothing)
        if length(idx) == (n - 1)
            df_BLUEs[idx, "Name"] = replace.(df_BLUEs[idx, "Name"], "entries: " => "")
            df_BLUEs[1, "Name"] = setdiff(entries, df_BLUEs[idx, "Name"])[1]
            df_BLUEs[2:end, "Coef."] .+= df_BLUEs[1, "Coef."]
            DataFrame(entries = entries)
        end
        # BLUPs
        df_BLUPs = DataFrame(only(raneftables(model)))
        # Find the columns referring to the nested entry effects, and excluding the spatial residual effects
        idx_col = findall(match.(r"blocks|rows|cols", names(df_BLUPs)) .== nothing)
        df_BLUPs = df_BLUPs[:, idx_col]

        additive_genetic = simulated_effects[1].additive_genetic
        for i = 2:length(simulated_effects)
            additive_genetic += simulated_effects[i].additive_genetic
        end
        additive_genetic ./= length(simulated_effects)
        cor(hcat(sum.(eachrow(df_BLUPs[:, 2:end])), additive_genetic))


    end

    false
end
