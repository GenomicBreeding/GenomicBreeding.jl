using Pkg
using JuliaFormatter
Pkg.update()
Pkg.activate(".")
# Format
format(".")
# Test
include("runtests.jl")
