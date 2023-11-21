using LinearAlgebra
using SparseArrays
using Plots
plotly()

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
    stored_errors = zeros(Float64,K);
    x = spzeros(Float64,length(A[:,1]))
    norm_b = norm(b)
    t = time()
    for i in 1:K
        x .= M\(T*x+b);
        #x .= diag(M).^(-1).*(T*x+b)
        rel_err = norm(A*x-b)/norm_b
        stored_errors[i] = rel_err
        if rel_err <= TOL
            println("Finished at $i iterations")
            break
        end
        if i == K
            println("Did not converge before maxit")
        end
    end
    dt1 = time()-t;

    println("Computation time Jacobi: $dt1 seconds")

    return x, stored_errors[stored_errors .!=0]
end

#Conjugate method.
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
        Ap = A*p;
        a = dot(p, r) / dot(p, Ap)
        x .= x .+ a .* p
        r_old .= r
        r .= r .- a .* (Ap)

        rel_err = norm(A * x - b) / norm_b
        stored_errors[i+1] = rel_err
        if  rel_err < TOL
            println("Finished at $i iterations")
            break
        end
    end
    dt = time() - t;
    println("Computation time cgm: $dt seconds")
    return x, stored_errors[stored_errors .!= 0]
end


function analysis(method, K, n_list, d_list, TOL = 1e-8, linearcomp = false)

    # method:: Int: 1=Jacobi, 2=CG
    # K:: Int: max_iterations
    # n:: Vector: discretization resolution
    # d:: Vector: dimension
    
    # For clear example where iterative methods are superior: n = 15, d = 4
    # Also, when excluding Jacobi, one can run n = 40, d = 3
    plot1 = plot()
    for n in n_list
        for d in d_list
            N = n^d
            b = sprand(N,1.0,rand,Float64)
            A = lap(n,d)
    
            if method == 1
                println("\nn = $n, d = $d. Computing using Jacobi...")
                x, SE = jacobi(A,b,K,TOL);
            end
            if method == 2
                println("\nn = $n, d = $d. Computing using Conjugate gradient...")
                x, SE = cgm(A, b, K, TOL);
            end
    
            err = norm(A*x-b)/norm(b)
            println("Relative error: $err")
    
            plot!(SE, yscale=:log10, label = "n = $n, d = $d")
            xlabel!("iteration")
            ylabel!("L_2 norm, relative error")
    
            if linearcomp == true
                println("\nn = $n, d = $d. Computing using Linear Solver...")
                t = time()
                A\Vector(b);
                dt2 = time()-t
                println("Computation time LS: $dt2 seconds")
            end
    
        end
    end
    display(plot1)
    end