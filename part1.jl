using LinearAlgebra
using SparseArrays
using MAT
using Plots
include("lap.jl")

#t = time() dt = time() - t

"""
mat = matopen("SF2520/SF2520-Lab5/convdiff.mat")
A = read(mat)["A"]
close(mat)
"""

n = 100
d = 2
N = n^d

b = sprand(N,1.0,rand,Float32)

K = 10000

A = lap(n,d)

M = spdiagm(diag(A)); T = M-A;

x = spzeros(N)
stored_errors = zeros(K);
TOL = 0.01
t = time()
for i in 1:K
    global x
    global rel_err
    x = M\(T*x+b);
    rel_err = norm(A*x-b)/norm(b)
    stored_errors[i] = rel_err
    if rel_err <= TOL
        println("Finished at $i iterations")
        break
    end
end
dt1 = time()-t;

println("Iterative time: $dt1 seconds, relative error: $rel_err")
plot(1:K, stored_errors, yscale=:log10)


t = time()
x1 = A\Vector(b)
dt2 = time()-t


println("Linear solver time: $dt2 seconds")
