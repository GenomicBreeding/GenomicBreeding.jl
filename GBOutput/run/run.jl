i = ARGS[1]
using GenomicBreeding
import GenomicBreeding: ols, ridge, lasso, bayesa, bayesb, bayesc
input = readjld2(GBInput, fname = joinpath("GBOutput/run", string("GBInput-", i, ".jld2")))
output = analysis(input)
display(output)
