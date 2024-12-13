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

# TODO: mixed model analysis
# NOTE@: 
# (1) Avoid over-parameterisation we'll have enough of that with the genomic prediction models
# (2) We will fit mixed models
# (3) We prefer REML over ML
# (4) We will compare BLUEs vs BLUPs of entries
function analyse(trials::Trials)::Bool
    # trials, _ = simulatetrials(genomes = simulategenomes());
    # trials, _ = simulatetrials(genomes = simulategenomes(n=5), n_years=2, n_seasons=2, n_harvests=1, n_sites=2, n_replications=10); trials.rows .= "Row_1"; trials.cols .= "Col_1";
    df::DataFrame = tabularise(trials)
    explanatory_variables = [
        "years",
        "seasons",
        "harvests",
        "sites",
        "replications",
        "blocks",
        "rows",
        "cols",
        "entries",
        "populations",
    ]
    explanatory_variables = filter(x -> length(unique(df[!, x])) > 1, explanatory_variables)
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
    # Total number of factorial combinations
    n_combin::Int64 = 1
    for v in explanatory_variables
        n_combin *= length(unique(df[!, v]))
    end

    trait = unique(trials.traits)[1]

    if (n_combin <= 500) && (n_combin <= nrow(df))
        # Includes saturated models - not too much sources of variation and dataset is tall instead of wide
        # f = @eval(@formula($(Meta.parse(string(traits[i], " ~ 1 + ", join(explanatory_variables, "*"))))))
        # # OLS
        # X = modelmatrix(f, df)
        # y = df[!, traits[i]]
        # _, sv = coefnames(apply_schema(f, schema(f, df)))
        # V = try 
        #     inv(X'*X)
        # catch
        #     pinv(X'*X)
        # end
        # β_ols = V*X'*y
        # ϵ_ols = y - X*β_ols
        # σ²_ϵ_ols = var(ϵ_ols)
        # σ²_β_ols = σ²_ϵ_ols .* diag(V)
        # df_ols = DataFrame(sv=sv, β_ols=β_ols, σ²_β_ols=σ²_β_ols)
        nester_vars::Array{String,1} =
            setdiff(explanatory_variables, ["entries", "blocks", "rows", "cols"])
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

        vec_formulae = [
            @eval(
                @formula(
                    $(Meta.parse(string(trait, " ~ 1 + years + (1 + entries|years)")))
                )
            ),
            @eval(
                @formula(
                    $(Meta.parse(
                        string(trait, " ~ 1 + years + zerocorr(1 + entries|years)"),
                    ))
                )
            ),
        ]

        f = @eval(
            @formula(
                $(Meta.parse(
                    string(trait, " ~ 1 + years + (1 + entries|(seasons:years))"),
                ))
            )
        )

        m_REML = fit(MixedModel, f, df, REML = true)
        deviance(m_REML)
        VarCorr(m_REML)
        df_BLUPs = DataFrame(only(raneftables(m_REML)))

    else
        # Does not include saturated models - too much sources of variation
        f = @eval(
            @formula(
                $(Meta.parse(string(trait, " ~ 1 + ", join(explanatory_variables, "*"))))
            )
        )
        X = modelmatrix(f, df)
    end

    false
end
