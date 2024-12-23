function performance(y_true::Vector{Float64}, y_pred::Vector{Float64})::Dict{String, Float64}
    # y_true = rand(10)
    # y_pred = rand(10)
    if var(y_true) < 1e-12
        @warn "No variance in y_true."
    end
    if var(y_pred) < 1e-12
        @warn "No variance in y_pred."
    end

    UnicodePlots.scatterplot(y_true, y_pred)

    diff::Vector{Float64} = y_true - y_pred
    mbe::Float64 = mean(diff)
    mae::Float64 = mean(abs.(diff))
    mse::Float64 = mean(diff.^2)
    rmse::Float64 = sqrt(mse)
    ρ_pearson::Float64 = cor(y_true, y_pred)
    ρ_spearman::Float64 = StatsBase.corspearman(y_true, y_pred)
    ρ_kendall::Float64 = StatsBase.corkendall(y_true, y_pred)
    r²::Float64 = 1.00 - sum(diff.^2) / sum((y_true .- mean(y_true)).^2)
    Dict(
        "mbe" => mbe,
        "mae" => mae,
        "mse" => mse,
        "rmse" => rmse,
        "ρ_pearson" => ρ_pearson,
        "ρ_spearman" => ρ_spearman,
        "ρ_kendall" => ρ_kendall,
        "r²" => r²,
    )
end

