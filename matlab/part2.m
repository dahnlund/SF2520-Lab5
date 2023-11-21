%% Computer Exercise 5
clear;clc

mat = load('../data/cooling_flange.mat');
A = mat.("A");

N = length(A(:,1));
b = rand(N,1);

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

%% pcg, M as ichog

L = ichol(A);
tic()
[x_pcgc,FLAG,RELRES,ITER,RESVEC_pcgc] = pcg(A, b, 1e-04, 10000, L, L');
time_pcgc = toc();
fprintf("pcg M=diag::      Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.0f, Flag: %.0f\n", ...
    time_pcgc, RELRES, ITER, FLAG)


%% Compare convergence
semilogy(RESVEC_cg);hold on; semilogy(RESVEC_pcgd); semilogy(RESVEC_pcgc)
xlabel("Iterations")
ylabel("RELRES")
legend("PCG", "PCG (diag)", "PCG (ichol)")

%% Part C

mat = load('../data/convdiff.mat');
A = mat.("A");

N = length(A(:,1));
b = rand(N,1);

mat = load('../data/b.mat');
b = mat.("b")';
mat = load('../data/x_gm.mat');
x_gm = mat.("x_gm")';
disp(norm(A*x_gm-b)/norm(b))

M = ilu(A);
tic()
[x_gm,FLAG,RELRES,ITER,RESVEC_gm] = gmres(A, b, false, 1e-4, 10000, M);
time_gm = toc();
fprintf("\nGMRES::     Computation time: %.04f seconds. RELRES: %.05d, Iterations: %.0f, Flag: %.0f\n",time_gm, RELRES, ITER, FLAG)





