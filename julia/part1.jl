using LinearAlgebra
using SparseArrays
#using MAT
using Plots
using IterativeSolvers
using IncompleteLU

include("src.jl")

#t = time() dt = time() - t

"""
mat = matopen("SF2520/SF2520-Lab5/convdiff.mat")
A = read(mat)["A"]
close(mat)
"""

n = 10
d = 4
N = n^d
K = 1000
TOL = 1e-8

b = sprand(N,1.0,rand,Float64)

A = lap(n,d)

x, SE1 = jacobi(A,b,K,TOL);


t = time()
x1 = A\Vector(b)
dt2 = time()-t
println("Linear solver time: $dt2 seconds")

TOL = 1e-8
x_conj, SE2 = cgm(A, b, K, TOL);

err_j = norm(A*x-b)/norm(b)
println("\nRelative error Jacobi: $err_j")
err_cgm = norm(A*x_conj-b)/norm(b)
println("\nRelative error conjugate descent: $err_cgm")
err_ls = norm(A*x1-b)/norm(b)
println("\nRelative error Linear Solver: $err_ls")

plot1 = plot(SE1, label="Jacobi convergence", yscale=:log10)
plot!(SE2, label="Conjugate descent convergence", yscale=:log10)
xlabel!("iteration")
ylabel!("L_2 norm, relative error")
display(plot1)

P = ilu(A, τ = 0.1);
x3 = bicgstabl(A,b,2, Pl = P, max_mv_products = 2000)
x4 = cg(A,b)


println(norm(A*x3-b)/norm(b))
println(norm(A*x4-b)/norm(b))