using LinearAlgebra
using SparseArrays
using Plots
include("src.jl")




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
