
n =50;
d = 3;
N = n^d;

b = rand(N,1);

K = 1000;

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

semilogy(1:K, stored_errors)


tic();
x1 = A\b;
dt2 = toc();

%fprintf("Iterative time: $dt1 seconds, TOL = $TOL")
%fprintf("Linear solver time: $dt2 seconds")

