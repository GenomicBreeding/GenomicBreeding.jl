function partitionkfolds(;genomes::Genomes, phenomes::Phenomes, k::Int64=10, r::Int64=3)::Vector{CV} end

function partitionpairwisepops(;genomes::Genomes, phenomes::Phenomes, r::Int64=1)::CV end

function partitionleaveoneoutpops(;genomes::Genomes, phenomes::Phenomes, r::Int64=1)::CV end

function cv!(cvs::Vector{CV}; genomes::Genomes, phenomes::Phenomes)::Nothing end