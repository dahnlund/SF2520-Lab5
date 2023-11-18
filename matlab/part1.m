
n =40;
d = 2;
N = n^d;

b = rand(N,1);

K = 10000;

A = lap(n,d);

M = spdiags(diag(A),0,N,N); T = M-A;

x = zeros(N,1);
stored_errors = zeros(K);
TOL = 0.01;
tic();
for i = 1:K
    x = M\(T*x+b);
    rel_err = norm(A*x-b)/norm(b);
    stored_errors(i) = rel_err;
    if rel_err <= TOL
        disp("Finished at $i iterations")
        break
    end
end
dt1 = toc();
fprintf("Iterative time: %.04f seconds, relative error: %.04f: \n", dt1,rel_err)
semilogy(1:i, stored_errors(1:i))


tic();
x1 = A\b;
dt2 = toc();

fprintf("Linear solver time: %.04f seconds\n", dt2)

%% Conjugate method

n =100;
d = 3;
N = n^d;

b = rand(N,1);

x = zeros(N,1);
K = 1000;

A = lap(n,d);

M = spdiags(diag(A),0,N,N); T = M-A;


p = x;
r = b;
TOL = 1e-8;
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
    if norm(A*x-b)/norm(b) < TOL
        break
    end
end

