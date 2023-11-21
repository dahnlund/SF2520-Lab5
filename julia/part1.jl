using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


K = 1000;
n = [10 12 15];
d = [2 3 4];
TOL = 1e-6
linearcomp = true

analysis(1, K, n, d, TOL)
analysis(2, K, n, d, TOL, linearcomp) # Run analysis for CG