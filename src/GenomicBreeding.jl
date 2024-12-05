module GenomicBreeding

include("core/genomes.jl")
include("core/phenomes.jl")
include("core/trials.jl")
include("simulation/simulate_genomes.jl")
include("simulation/simulate_trials.jl")

export Genomes, Phenomes, Trials
export checkdims, dimensions
export simulategenomes, simulategenomiceffects, simulatetrials

end
