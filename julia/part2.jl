using MAT
using LinearAlgebra
using SparseArrays
using Plots
using IterativeSolvers
using IncompleteLU


mat = matopen("data/cooling_flange.mat")
A = read(mat)["A"]
close(mat)
N = length(A[:,1])
"""
P = ilu(A, τ = 0.1);
x3 = bicgstabl(A,b,2, Pl = P, max_mv_products = 2000)
x4 = cg(A,b)


println(norm(A*x3-b)/norm(b))
println(norm(A*x4-b)/norm(b))

"""

b = sprand(N,1.0,rand,Float64)

t = time()
x_ls = A\Vector(b);
time_ls = time()-t

err_ls = norm(A*x_ls-b)/norm(b)
println("Computation time Linear Solver: $time_ls. RELRES: $err_ls")

"""
t = time()
x_cg, RESVEC = cg(A, b, reltol = 1e-1 , log = true)
time_cg = time()-t

err_cg = norm(A*x_cg-b)/norm(b)
println("Computation time Linear Solver: $time_cg. RELRES: $err_cg")
"""


LU = ilu(A,τ = 0.1)
t = time()
x_pcg = cg(A, b, Pl = LU, reltol = 1e-4)
time_pcg = time()-t

err_pcg = norm(A*x_pcg-b)/norm(b)
println("Computation time Linear Solver: $time_pcg. RELRES: $err_pcg")
