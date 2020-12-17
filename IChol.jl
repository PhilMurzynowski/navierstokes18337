using Printf
using Kronecker
using LinearAlgebra

"""
Incomplete Cholesky factorization
Do not need to optimize as only need to get a sense of how it looks
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


function genPoissonStreamMtx(size, opts)
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
    a = genPoissonStreamMtx(opts["N"], opts)
    @printf "a\n"
    display(a)
    @printf "L\n"
    L = ichol(a)
    display(L)
    @printf "Linv"
    Linv = inv(L)
    display(Linv)
    @printf "Minv\n"
    Minv = Linv'*Linv
    display(Minv)
end

#N = 3
#opts = Dict("N"=>N,
#            "h"=>1/N
#            )
#tests(opts)