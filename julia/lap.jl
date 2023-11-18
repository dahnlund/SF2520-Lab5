using LinearAlgebra
using SparseArrays

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


