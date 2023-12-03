
n =12;
d = 4;
N = n^d;

b = rand(N,1);

K = 10000;

A = lap(n,d);
bnorm = norm(b);
M = spdiags(diag(A),0,N,N); T = M-A;

x = zeros(N,1);
stored_errors = zeros(K,1);
TOL = 1e-10;
tic();
for i = 1:K
    x = M\(T*x+b);
    rel_err = norm(A*x-b)/bnorm;
    stored_errors(i) = rel_err;
    if rel_err <= TOL
        fprintf("Finished at %.0f iterations\n", i)
        break
    end
end
dt1 = toc();
fprintf("Iterative time: %.04f seconds, relative error: %.04d: \n", dt1,rel_err)
semilogy(1:i, stored_errors(1:i))


tic();
x1 = A\b;
dt2 = toc();

fprintf("Linear solver time: %.04f seconds\n", dt2)

%% Conjugate method

n =750;
d = 1;
N = n^d;

b = rand(N,1);

x = zeros(N,1);
K = 10000;

A = lap(n,d);
norm_b = norm(b);
M = spdiags(diag(A),0,N,N); T = M-A;

stored_errors = zeros(K,1);
p = x;
r = b;
TOL = 1e-10;
tic();
for i = 0:K
    if i == 0
        beta = 0;
    else
        beta = r'*r/(r_old'*r_old);
    end
    p = r + beta*p;
    a = p'*r/(p'*A*p);
    
    x = x+a*p;
    r_old = r;
    r = r - a*A*p;
    rel_err = norm(A*x-b)/norm_b;
    stored_errors(i+1) = rel_err;
    if rel_err < TOL
        fprintf("Finished at %.0f iterations\n", i)
        break
    end
end
dt3 = toc();
fprintf("Iterative time CG: %.04f seconds, relative error: %.04d: \n", dt3,rel_err)
semilogy(1:i, stored_errors(1:i))