using LinearAlgebra
using SparseArrays
using MAT
using Plots
include("src.jl")

#t = time() dt = time() - t

"""
mat = matopen("SF2520/SF2520-Lab5/convdiff.mat")
A = read(mat)["A"]
close(mat)
"""

n = 60
d = 3
N = n^d
K = 1000
TOL = 1e-5

b = sprand(N,1.0,rand,Float64)

A = lap(n,d)

x, SE1 = jacobi(A,b,K,TOL);


t = time()
x1 = A\Vector(b)
dt2 = time()-t
println("Linear solver time: $dt2 seconds")

TOL = 1e-5
x_conj, SE2 = cgm(A, b, K, TOL);

err_j = norm(A*x-b)/norm(b)
println("\nRelative error Jacobi: $err_j")
err_cgm = norm(A*x_conj-b)/norm(b)
println("\nRelative error conjugate descent: $err_cgm")
err_ls = norm(A*x1-b)/norm(b)
println("\nRelative error Linear Solver: $err_ls")