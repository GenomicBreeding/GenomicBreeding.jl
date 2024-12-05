using Pkg
using JuliaFormatter
Pkg.activate(".")
format(".")

include("runtests.jl")
