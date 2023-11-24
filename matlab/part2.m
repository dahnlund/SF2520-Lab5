%% Computer Exercise 5
clear;clc

mat = load('../data/cooling_flange.mat');
A = mat.("A");

N = length(A(:,1));
b = rand(N,1);
norm_b = norm(b);

%% Linear Solver
tic
x_ls = A\b;
time_ls = toc();

err_ls = norm(A*x_ls-b)/norm(b);
fprintf("Linear Solver::      Computation time: %.04f seconds. RELRES: %.05d\n", time_ls, err_ls)

%% Regular pcg
tic()
[x_cg,FLAG,RELRES,ITER,RESVEC_cg] = pcg(A, b, 1e-04, 10000);
time_cg = toc();
fprintf("Regular pcg::      Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.0f, Flag: %.0f\n", time_cg, RELRES, ITER, FLAG)

%% pcg, M as diagnonal

M = spdiags(diag(A),0,N,N);
tic()
[x_pcgd,FLAG,RELRES,ITER,RESVEC_pcgd] = pcg(A, b, 1e-04, 10000, M);
time_pcgd = toc();
fprintf("pcg M=diag::      Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.0f, Flag: %.0f\n", ...
    time_pcgd, RELRES, ITER, FLAG)

%% pcg, M as ichol

L = ichol(A);
tic()
[x_pcgc,FLAG,RELRES,ITER,RESVEC_pcgc] = pcg(A, b, 1e-04, 10000, L, L');
time_pcgc = toc();
fprintf("pcg M=ichol::      Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.0f, Flag: %.0f\n", ...
    time_pcgc, RELRES, ITER, FLAG)


%% Compare convergence
semilogy(RESVEC_cg/norm_b);hold on; semilogy(RESVEC_pcgd/norm_b); semilogy(RESVEC_pcgc/norm_b)
xlabel("Iterations")
ylabel("RELRES")
legend("PCG", "PCG (diag)", "PCG (ichol)")

%% Part C

mat = load('../data/convdiff.mat');
A = mat.("A");

N = length(A(:,1));
b = rand(N,1);
norm_b = norm(b);

%% Solve with GE
tic
x_gmls = A\b;
time_gmls = toc();

err_gmls = norm(A*x_gmls-b)/norm(b);
fprintf("Linear Solver::      Computation time: %.04f seconds. RELRES: %.05d\n", time_gmls, err_gmls)

%% Using GMRES

[L,U] = ilu(A);
tic()
[x_gm,FLAG,RELRES,ITER,RESVEC_gm] = gmres(A, b, N, 1e-5, 3000, L,U);
time_gm = toc();
fprintf("GMRES::  Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.01f, Flag: %.01f\n",time_gm, RELRES, ITER(2), FLAG)
semilogy(RESVEC_gm/norm_b); hold on

tic()
[x_gmw,FLAGw,RELRESw,ITERw,RESVEC_gmw] = gmres(A, b, N, 1e-5, 3000);
time_gmw = toc();
fprintf("GMRES no precond.::  Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.01f, Flag: %.01f\n",time_gmw, RELRESw, ITERw(2), FLAGw)
semilogy(RESVEC_gmw/norm_b)



