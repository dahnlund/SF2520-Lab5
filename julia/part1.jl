using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")


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
