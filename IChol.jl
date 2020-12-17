using Printf
using Kronecker
using LinearAlgebra

"""
Incomplete Cholesky factorization
Not currently optimized as only need to ever call it once.
Would be very good to specify it as a LowerTriangular type
Disclaimer: Not my code
Source: wikipedia => https://en.wikipedia.org/wiki/Incomplete_Cholesky_factorization
"""
function ichol(a)
	n = size(a)[1];

	for k=1:n
		a[k,k] = sqrt(a[k,k]);
		for i in k+1:n
		    if a[i,k]!=0
		        a[i,k] = a[i,k]/a[k,k]           
            end
		end
		for j in k+1:n
            for i in j:n
                # check if nonzero
		        if a[i,j]!=0
		            a[i,j] = a[i,j]-a[i,k]*a[j,k] 
		        end
		    end
		end
	end

    # zero out upper triangle
    for i=1:n
        for j=i+1:n
            a[i,j] = 0
        end
    end            
    return a
end


function genPoissonMtx(size, opts)
    N = size
    hi = 1/opts["h"]

    L = zeros(N, N)
    L[diagind(L, 0)] .= 2*hi^2
    L[diagind(L, -1)] .= -1*hi^2
    L[diagind(L, 1)] .= -1*hi^2
    # build Poisson pressure matrix with Kronecker
    # cant use I with Kronecker, not working
    Id = 1*Matrix(I, N, N)
    P = (Id ⊗ L) + (L ⊗ Id)
    BC = zeros(N, N)
    return P
end

function tests(opts)
    @printf "test start\n\n\n"
    A = genPoissonMtx(opts["N"], opts)
    @printf "A\n"
    display(A)
    @printf "L\n"
    L = ichol(A)
    display(L)
    @printf "Linv"
    Linv = inv(L)
    display(Linv)
    #@printf "Minv\n"
    #Minv = Linv'*Linv
    #display(Minv)
    @printf "LinvALinv'\n"
    A_new = Linv*A*Linv'
    display(A_new)
end

N = 5
opts = Dict("N"=>N,
            "h"=>1/N
            )
tests(opts)