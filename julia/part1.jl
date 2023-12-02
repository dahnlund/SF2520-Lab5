using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


K = 30000;
n = [10^4];
d = [1];

comb = false
TOL = 1e-10
linearcomp = true

plot1 = analysis(1, K, n, d, comb, TOL)
display(plot1)

#plot2 = analysis(2, K, n, d, comb, TOL, linearcomp, plot1) # Run analysis for CG
#display(plot2)