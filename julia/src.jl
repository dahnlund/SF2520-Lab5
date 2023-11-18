using LinearAlgebra
using SparseArrays
using Plots

function lap(n,d)
    # LAP
    #    A = lap(N,D) returns the system matrix corresponding to the 
    #    standard second order discretization of the Laplace operator 
    #    for D dimensions and N unknowns in each coordinate direction.
    #   The size of the matrix is thus N^D times N^D.
    
    e = ones(n);
    A1 = -spdiagm(-1 => e[1:end-1], 0 => (-2*e), 1=> e[1:end-1])
    I1 = sparse(I,n,n)
    
    A=A1;
    I_=I1;
    for _ in 2:d
      A = kron(A,I1)+kron(I_,A1);
      I_ = kron(I_,I1);
    end
    A = A*n^2;
    return A
end


#Jacobi Method
function jacobi(A, b, K, TOL)
    M = spdiagm(diag(A)); T = M-A;
    stored_errors = zeros(K);
    x = spzeros(length(A[:,1]))
    norm_b = norm(b)
    t = time()
    for i in 1:K
        x .= M\(T*x+b);
        rel_err = norm(A*x-b)/norm_b
        stored_errors[i] = rel_err
        if rel_err <= TOL
            println("Finished at $i iterations")
            break
        end
    end
    dt1 = time()-t;

    println("Computation time Jacobi: $dt1 seconds")
    return x, stored_errors
end

#Conjugate method
function cgm(A, b, K, TOL)
    x = spzeros(length(A[:,1]))
    norm_b = norm(b)
    p = copy(x)
    stored_errors = zeros(K);
    r = copy(b)
    r_old = copy(r)
    t = time()
    for i in 0:K-1
        if i == 0
            beta = 0
        else
            beta = dot(r, r) / dot(r_old, r_old)
        end

        p .= r .+ beta .* p
        a = dot(p, r) / dot(p, A * p)
        x .= x .+ a .* p
        r_old .= r
        r .= r .- a .* (A * p)

        rel_err = norm(A * x - b) / norm_b
        stored_errors[i+1] = rel_err
        if  rel_err < TOL
            println("Finished at $i iterations")
            break
        end
    end
    dt = time() - t;
    println("Computation time cgm: $dt seconds")
    return x, stored_errors
end