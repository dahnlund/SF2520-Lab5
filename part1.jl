using LinearAlgebra
using SparseArrays
using MAT
using Plots
include("lap.jl")

#t = time() dt = time() - t

"""
mat = matopen("SF2520/SF2520-Lab5/convdiff.mat")
A = read(mat)["A"]
close(mat)
"""

n = 10
d = 1

b = rand(n,1)

K = 100

A = lap(n,d)

M = spdiagm(diag(A)); T = M-A;

x = spzeros(n)
for _ in 1:K
    global x
    x = M\(T*x+b);
end

println(norm(x-A\b))