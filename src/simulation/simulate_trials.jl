using Random
using StatsBase
using Distributions
using LinearAlgebra
using ProgressMeter
using Dates
using Plots
using StatsPlots
# Plots.backend(:plotly)
Plots.backend(:gr)

"""
# Simulate genomic effects

Simulate additive, dominance, and epistatic effects

## Arguments
- `genomes`: Genome struct includes the `n` entries x `p` loci-alleles combinations (`p` = `l` loci x `a-1` alleles)
- `f_additive`: proportion of the `l` loci with non-zero additive effects on the phenotype
- `f_dominance`: proportion of the `l*f_additive` additive effects loci with additional dominance effects
- `f_epistasis`: proportion of the `l*f_additive` additive effects loci with additional dominance effects

## Outputs

## Examples
# ```jldoctest; setup = :(using GenomicBreeding)
# ```
"""
function simulategenomiceffects(;
    genomes::Genomes,
    f_additive::Real = 0.01,
    f_dominance::Real = 0.10,
    f_epistasis::Real = 0.05,
    seed::Int = 42,
)::Tuple{Array{Real,1},Array{Real,1},Array{Real,1}}
    # genomes::Genomes = simulategenomes(n=100, l=2_000, n_alleles=3); f_additive::Real = 0.01; f_dominance::Real = 0.10; f_epistasis::Real = 0.05; seed::Int = 42;
    # Argument checks
    if !checkdims(genomes)
        throw(ArgumentError("simulategenomiceffects: error in the genomes input"))
    end
    if (f_additive < 0.0) || (f_additive > 1.0)
        throw(ArgumentError("We accept `f_additive` from 0.00 to 1.00."))
    end
    if (f_dominance < 0.0) || (f_dominance > 1.0)
        throw(ArgumentError("We accept `f_dominance` from 0.00 to 1.00."))
    end
    if (f_epistasis < 0.0) || (f_epistasis > 1.0)
        throw(ArgumentError("We accept `f_epistasis` from 0.00 to 1.00."))
    end
    # Genomes dimensions
    n::Int, p::Int, l::Int, max_n_alleles::Int = dimensions(genomes)
    # Number of loci with additive, dominance, and epistasis effects (minimum values of 1, 0, and 0 or 2 (if there is non-zero loci with epistasis then we expect to have at least 2 loci interacting), respectively; Note that these are also imposed in the above arguments checks)
    a::Int = Int(maximum([1, round(l * f_additive)]))
    d::Int = Int(maximum([0, round(l * f_additive * f_dominance)]))
    e::Int = Int(maximum([0, round(l * f_additive * f_epistasis)]))
    if e == 1
        e = 2 # if there is one epistatic locus then it should interact with at least one other locus
    end
    # Instantiate the output vectors
    β::Array{Real,1} = fill(0.0, p)
    δ::Array{Real,1} = fill(0.0, p)
    ξ::Array{Real,1} = fill(0.0, p)
    
    # Set randomisation seed
    rng::TaskLocalRNG = Random.seed!(seed)
    # Define the locations of the loci with non-zero genetic effects
    idx_additive::Array{Int,1} =
        StatsBase.sample(rng, 1:l, a, replace = false, ordered = true)
    idx_dominance::Array{Int,1} =
        StatsBase.sample(rng, idx_additive, d, replace = false, ordered = true)
    idx_epistasis::Array{Int,1} =
        StatsBase.sample(rng, idx_additive, e, replace = false, ordered = true)
    # Sample additive allele effects from a multivariate normal distribution
    μ_β_dist::Exponential = Distributions.Exponential(1.0)
    μ_β::Array{Real,1} = rand(rng, μ_β_dist, a)
    Σ_β::Array{Real,2} = μ_β * μ_β'
    for _ = 1:10
        if abs(det(Σ_β)) < 1e-12
            Σ_β[diagind(Σ_β)] .+= 1.0
        else
            break
        end
    end
    if abs(det(Σ_β)) < 1e-12
        throw(
            ArgumentError(
                "Please increase/reduce the proportion of loci with additive effects.",
            ),
        )
    end
    β_dist::MvNormal = Distributions.MvNormal(μ_β, Σ_β)
    B::Array{Real,2} = rand(rng, β_dist, max_n_alleles-1)
    # Define the loci-alleles combination indexes corresponding to the additive allele loci
    idx_p_additive = vcat((idx_additive*(max_n_alleles-1)) .- 1, idx_additive*(max_n_alleles-1))
    sort!(idx_p_additive)
    # Update the addivite allele effects
    β[idx_p_additive] = reshape(B', (a*(max_n_alleles-1), 1))

    




    (β, δ, ξ)

end


"""
# Simulate trials

## Arguments

## Output

## Examples
"""
function simulatetrials(;
    genomes::Genomes,
    n_traits::Int = 5,
    n_years::Int = 2,
    n_seasons::Int = 4,
    n_harvests::Int = 2,
    n_sites::Int = 4,
    n_replications::Int = 2,
    n_blocks::Int = 2,
    n_rows::Int = 2,
    n_cols::Int = 2,
)::Trials
    # phenotypes::Array{Union{Real,Missing},2}
    # traits::Array{String,1}
    # years::Array{String,1}
    # seasons::Array{String,1}
    # harvests::Array{String,1}
    # sites::Array{String,1}
    # replications::Array{String,1}
    # blocks::Array{String,1}
    # rows::Array{String,1}
    # cols::Array{String,1}
    # entries::Array{String,1}
    # populations::Array{String,1}
    Trials(n = 2, t = 3)
end
