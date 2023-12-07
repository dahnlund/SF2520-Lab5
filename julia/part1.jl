#CE5 SE2520. Authors: David Ahnlund, Emil Gestsson

using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


K = 10000;
n = [500 5000 30 10];
d = [1 1 3 4];

comb = false
TOL = 1e-4
linearcomp = true

plot1 = analysis(1, K, n, d, comb, TOL)
display(plot1)

plot2 = analysis(2, K, n, d, comb, TOL, linearcomp) # Run analysis for CG
display(plot2)