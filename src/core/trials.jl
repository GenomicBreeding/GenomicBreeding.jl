"""
# Trials struct 
    
Contains phenotype data across years, seasons, harvest, sites, populations, replications, blocks, rows, and columns

## Constructor

```julia
Trials(; n::Int64 = 2, p::Int64 = 2)
```

## Fields
- `phenotypes`: `n x t` matrix of numeric phenotype data which can have missing values
- `traits`: names of the traits `t` traits
- `years`: names of the years corresponding to each row in the phenotype matrix
- `seasons`: names of the seasons corresponding to each row in the phenotype matrix
- `harvests`: names of the harvests corresponding to each row in the phenotype matrix
- `sites`: names of the sites corresponding to each row in the phenotype matrix
- `replications`: names of the replications corresponding to each row in the phenotype matrix
- `blocks`: names of the blocks corresponding to each row in the phenotype matrix
- `rows`: names of the rows corresponding to each row in the phenotype matrix
- `cols`: names of the cols corresponding to each row in the phenotype matrix
- `entries`: names of the entries corresponding to each row in the phenotype matrix
- `populations`: names of the populations corresponding to each row in the phenotype matrix

## Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2)
Trials(Union{Missing, Float64}[missing missing], [#undef, #undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef], [#undef])

julia> fieldnames(Trials)
(:phenotypes, :traits, :years, :seasons, :harvests, :sites, :replications, :blocks, :rows, :cols, :entries, :populations)
```
"""
mutable struct Trials
    phenotypes::Array{Union{Float64,Missing},2}
    traits::Array{String,1}
    years::Array{String,1}
    seasons::Array{String,1}
    harvests::Array{String,1}
    sites::Array{String,1}
    replications::Array{String,1}
    blocks::Array{String,1}
    rows::Array{String,1}
    cols::Array{String,1}
    entries::Array{String,1}
    populations::Array{String,1}
    function Trials(; n::Int64 = 2, t::Int64 = 2)
        new(
            fill(missing, n, t),
            Array{String,1}(undef, t),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
            Array{String,1}(undef, n),
        )
    end
end

"""
    checkdims(trials::Trials)::Bool

Check dimension compatibility of the fields of the Trials struct

# Examples
```jldoctest; setup = :(using GenomicBreeding)
julia> trials = Trials(n=1, t=2);

julia> trials.entries = ["entry_1"];

julia> checkdims(trials)
true

julia> trials.entries = ["entering_2_entries", "instead_of_just_1"];

julia> checkdims(trials)
false
```
"""
function checkdims(trials::Trials)::Bool
    n, t = size(trials.phenotypes)
    if (t != length(trials.traits)) ||
       (n != length(trials.years)) ||
       (n != length(trials.seasons)) ||
       (n != length(trials.harvests)) ||
       (n != length(trials.sites)) ||
       (n != length(trials.replications)) ||
       (n != length(trials.blocks)) ||
       (n != length(trials.rows)) ||
       (n != length(trials.cols)) ||
       (n != length(trials.entries)) ||
       (n != length(trials.populations))
        return false
    end
    return true
end

function tabularise(trials::Trials)::DataFrame
    # trials::Trials, _ = simulatetrials(genomes = simulategenomes());
    df_ids::DataFrame = DataFrame(
        id = 1:length(trials.years),
        years = trials.years,
        seasons = trials.seasons,
        harvests = trials.harvests,
        sites = trials.sites,
        replications = trials.replications,
        blocks = trials.blocks,
        rows = trials.rows,
        cols = trials.cols,
        entries = trials.entries,
        populations = trials.populations,
    )
    df_phe::DataFrame = DataFrame(trials.phenotypes, :auto)
    rename!(df_phe, trials.traits)
    df_phe.id = 1:length(trials.years)
    df = innerjoin(df_ids, df_phe, on = :id)
    df
end

# # TODO: plot raw data
# function plot(trials::Trials)::Bool
#     # trials, _ = simulatetrials(genomes = simulategenomes())
#     df::DataFrame = tabularise(trials)
#     cor(trials.phenotypes)
#     plt = corrplot(trials.phenotypes)


#     false
# end

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
    df::DataFrame = tabularise(trials)
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
    # Exclude replications (residuals estimator) and population (entry-specific and can be used instead of entries) in the formula
    explanatory_variables = filter(x -> x != "replications", explanatory_variables)
    explanatory_variables = filter(x -> x != "populations", explanatory_variables)

    # Number of entries whose BLUEs or BLUPs we wish to extract
    entries::Vector{String} = sort(unique(df[!, "entries"]))
    n::Int64 = length(entries)

    # Iterate across each trait
    for trait in unique(trials.traits)
        # trait = unique(trials.traits)[1]
        # Find nester variables and permutate
        nester_vars::Array{String,1} = setdiff(explanatory_variables, ["entries"])
        # setdiff(explanatory_variables, ["entries", "blocks", "rows", "cols"])
        vars::Array{Array{String,1},1} = [[x] for x in nester_vars]
        for nv in nester_vars
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
            # var = vars[2]

            # Sort by decreasing nester
            spatials = var[[findall(var .== v)[1] for v in spatial_variables ∩ var]]
            nesters = var[[findall(var .== v)[1] for v in nester_variables ∩ var]]

            # Define the spatial errors
            spatial_error_column::String = join(spatials, "_x_")
            if spatial_error_column != ""
                try
                    eval(Meta.parse(string("df[!, \"", spatial_error_column, "\"]")))
                catch
                    eval(
                        Meta.parse(
                            string(
                                "df[!, \"",
                                spatial_error_column,
                                "\"] = string.(",
                                join(string.("df[!, \"", spatials, "\"]"), ", \"_x_\", "),
                                ")",
                            ),
                        ),
                    )
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
                        eval(Meta.parse(string("df[!, \"", interaction_column, "\"]")))
                    catch
                        eval(
                            Meta.parse(
                                string(
                                    "df[!, \"",
                                    interaction_column,
                                    "\"] = string.(",
                                    join(
                                        string.("df[!, \"", random_interaction, "\"]"),
                                        ", \"_x_\", ",
                                    ),
                                    ")",
                                ),
                            ),
                        )
                    end
                    if length(
                        unique(
                            eval(Meta.parse(string("df[!, \"", interaction_column, "\"]"))),
                        ),
                    ) > max_levels
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


        models::Vector{LinearMixedModel{Float64}} =
            Vector{LinearMixedModel{Float64}}(undef, length(formulae))
        deviances::Vector{Float64} = fill(Inf64, length(formulae))
        # @time Threads.@threads for i in eachindex(formulae)
        @time Threads.@threads for i = 1:10
            # i = 5
            f = @eval(@formula($(Meta.parse(formulae[i]))))
            model = MixedModel(f, df)
            model.optsum.maxtime = max_time_per_model
            @time fit!(model, REML = true)
            models[i] = model
            deviances[i] = deviance(model)
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
        additive_genetic = simulated_effects[1].additive_genetic
        for i = 2:length(simulated_effects)
            additive_genetic += simulated_effects[i].additive_genetic
        end
        additive_genetic ./= length(simulated_effects)
        cor(hcat(sum.(eachrow(df_BLUPs[:, 2:end])), additive_genetic))


    end

    false
end


try
    @warn "test"
catch
    println("HEY")
end
