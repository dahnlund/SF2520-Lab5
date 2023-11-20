using MAT
using LinearAlgebra
using SparseArrays
using Plots
using IterativeSolvers
using IncompleteLU

"Part 2c"

mat = matopen("data/convdiff.mat")
#mat = matopen("data/cooling_flange.mat")
A = read(mat)["A"]
close(mat)
N = length(A[:,1])


b = sprand(N,1.0,rand,Float64)

#Linear Solve
t = time()
x_ls = A\Vector(b);
time_ls = time()-t

err_ls = norm(A*x_ls-b)/norm(b)
println("Computation time Linear Solver: $time_ls. RELRES: $err_ls")

# GMRES
LU = ilu(A, τ = 0.1)
t = time()
x_gm = bicgstabl(A,b, Pl = LU)
time_gm = time()-t

err_gm = norm(A*x_gm-b)/norm(b)
println("Computation time BiCG: $time_gm. RELRES: $err_gm")
plot1 = plot(RESVEC)

