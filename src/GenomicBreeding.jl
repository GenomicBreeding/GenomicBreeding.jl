module GenomicBreeding
include("core/genomes.jl")
include("core/phenomes.jl")
include("core/trials.jl")
include("simulation/simulate_genomes.jl")
export Genomes, Phenomes, Trials
export check
export simulategenomes
end
