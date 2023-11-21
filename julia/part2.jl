using MAT
using LinearAlgebra
using SparseArrays
using Plots
using IterativeSolvers
using IncompleteLU
using Preconditioners

"Part 2"
#mat = matopen("data/convdiff.mat")
mat = matopen("data/cooling_flange.mat")
A = read(mat)["A"]
close(mat)
N = length(A[:,1])
b = sprand(N,1.0,rand,Float64)

TOL = 1e-4
K = 1000;
#Linear Solve
t = time()
x_ls = A\Vector(b);
time_ls = time()-t

err_ls = norm(A*x_ls-b)/norm(b)
println("Computation time Linear Solver: $time_ls. RELRES: $err_ls")


# Regular PCG
t = time()
x_regpcg, RESVEC_regpcg = minres(A,b, false, false, 1e-4, 1e-4, K, true, false);
time_ls = time()-t
err_ls = norm(A*x_regpcg-b)/norm(b)
println("Computation time regular pcg: $time_regpcg. RELRES: $err_regpcg")





"""
# GMRES
LU = ilu(A, τ = 0.1)
t = time()
x_gm, RESVEC = bicgstabl(A,b, Pl = LU, log = true, reltol = 1e-4)
time_gm = time()-t

err_gm = norm(A*x_gm-b)/norm(b)
println("Computation time BiCG: $time_gm. RELRES: $err_gm")
plot1 = plot(RESVEC)


t = time()
x_gm, RESVEC = bicgstabl(A,b, log = true, reltol = 1e-4)
time_gm = time()-t

err_gm = norm(A*x_gm-b)/norm(b)
println("Computation time BiCG (without ilu): $time_gm. RELRES: $err_gm")
plot!(RESVEC)
display(plot1)

"""