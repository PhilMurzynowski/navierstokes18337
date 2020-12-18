using Printf
using Kronecker
using LinearAlgebra
using BandedMatrices, BlockBandedMatrices
using SparseArrays

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


function genPoissonMtx(size, spacing)
    N = size
    hi = 1/spacing

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

"""
Uses BandedBlock matrices for less space and more optimized operations.
Followed: https://discourse.julialang.org/t/best-way-to-construct-2-d-laplacian-banded-matrix-bandedmatrices-sparse-or-blockbandedmatrices/29123
"""
function genPoissonMtxBanded(size, spacing)
    n = size
    h = spacing
    L1D = BandedMatrix(0 => Fill(2/h^2,n), 1 => Fill(-1/h^2,n-1), -1 => Fill(-1/h^2,n-1))
    D_xx = BandedBlockBandedMatrix(kron(L1D, Eye(n)))
    D_yy = BandedBlockBandedMatrix(kron(Eye(n), L1D))
    P = D_xx + D_yy
    return P
end

#P = genPoissonMtxBanded(4, 1/4)
#display(P)