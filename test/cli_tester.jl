using Pkg
using JuliaFormatter
Pkg.activate(".")
Pkg.update()
# Format
format(".")
# Test
include("runtests.jl")
