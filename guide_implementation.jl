using Kronecker
using LinearAlgebra

N = 10
dx, dy = 1/N, 1/N
# if dx = dy, can abstract to h
h = dx
dxi, dyi, hi = 1/dx, 1/dy, 1/h
# NOT entirely original implementation here, testing guide
# what about boundary conditions, sentinel?
# the plus one might be handling? unclear

# interleave u and v update to mitigate cache misses for large matrices?

for j in jmin:jmax
    for i in imin+1:imax
        v_here = 0.25 * ( v[i-1, j] + v[i-1, j+1] + v[i, j] + v[i, j+1])
        update = (u[i-1,j] - 2*u[i, j] + u[i+1, j]) * dxi^2
        update += (u[i, j-1] - 2*u[i, j] + u[i, j+1]) * dyi^2
        update -= u[i, j] * (u[i+1, j] - u[i-1, j]) * 0.5 * dxi
        update -= v_here * (u[i, j+1] - u[i, j-1]) * 0.5 * dyi
        update *= dt
        u_new[i, j] += u[i, j] + update
    end
end

for j in jmin+1:jmax
    for i in imin:imax
        u_here = 0.25 * ( u[i, j-1] + u[i, j] + u[i+1, j-1] + u[i+1, j])
        update = (v[i-1,j] - 2*v[i, j] + v[i+1, j]) * dxi^2
        update += (v[i, j-1] - 2*v[i, j] + v[i, j+1]) * dyi^2
        update -= u_here * (v[i+1, j] - v[i-1, j]) * 0.5 * dxi
        update -= v[i, j] * (v[i, j+1] - u[i, j-1]) * 0.5 * dyi
        update *= dt
        v_new[i, j] += v[i, j] + update
    end
end

# can describe the general structure with Kronecker products
# with some small modifications due to Von Neumman boundary conditions
# and also first cell is kept constant as pressure is defined up to a constant
#   other cells in reference to cell that is kept constant
# L is tridiagonal Laplacian Matrix
L = zeros(N, N)
L[diagind(L, 0)] .= 2*hi^2
L[diagind(L, -1)] .= -1*hi^2
L[diagind(L, 1)] .= -1*hi^2
# build Poisson pressure matrix with Kronecker
# cant use I with Kronecker, not working
Id = 1*Matrix(I, N, N)
P = (Id ⊗ L) + (L ⊗ Id)
# Von Neumann BC and reference point modifications