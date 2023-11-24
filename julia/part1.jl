using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


K = 1000;
n = [10 12 15];
d = [2 3 4];

TOL = 1e-10
linearcomp = true

plot1 = analysis(1, K, n, d, TOL)
plot2 = analysis(2, K, n, d, TOL, linearcomp, plot1) # Run analysis for CG

display(plot1)
display(plot2)