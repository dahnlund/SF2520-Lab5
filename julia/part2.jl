using MAT
using LinearAlgebra
using SparseArrays
using Plots
using IterativeSolvers
using IncompleteLU
using Preconditioners
using LinearSolve

"Part 2"
mat = matopen("data/convdiff.mat")
#mat = matopen("data/cooling_flange.mat")
A = read(mat)["A"]
close(mat)
N = length(A[:,1])
b = sprand(N,1.0,rand,Float64)

TOL = 1e-4
K = 1000;
#Linear Solve
#t = time()
#x_ls = A\Vector(b);
#time_ls = time()-t

#err_ls = norm(A*x_ls-b)/norm(b)
#println("Computation time Linear Solver: $time_ls. RELRES: $err_ls")

Pl  = CholeskyPreconditioner(A,5)
t = time()
prob = LinearProblem(A,b)
x_lsx = solve(prob, KrylovJL_CG(), Pl= Pl).u
time_lsx = time()-t
err_lsx = norm(A*x_lsx - b)/norm(b)
println("Computation time LinearSolve: $time_lsx. RELRES: $err_lsx")
#KrylovJL_BICGSTAB() best so far



t = time()
prob = LinearProblem(A,b)
x_lsx = solve(prob, KrylovJL_BICGSTAB()).u
time_lsx = time()-t
err_lsx = norm(A*x_lsx - b)/norm(b)
println("Computation time LinearSolve: $time_lsx. RELRES: $err_lsx")