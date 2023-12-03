using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


K = 10000;
n = [3000 4000 15 7];
d = [1 1 3 4];

comb = false
TOL = 1e-10
linearcomp = true

plot1 = analysis(1, K, n, d, comb, TOL, linearcomp)
display(plot1)

plot2 = analysis(2, K, n, d, comb, TOL) # Run analysis for CG
display(plot2)