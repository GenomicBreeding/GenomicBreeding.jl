using Pkg
using JuliaFormatter
Pkg.activate(".")
# Format
format(".")
# Test
include("runtests.jl")
