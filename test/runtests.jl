using Test
using Random
using LinearAlgebra

using DAESolve

include("pendulum.jl")

include("test_massmatrix.jl")
include("test_solvedae.jl")
