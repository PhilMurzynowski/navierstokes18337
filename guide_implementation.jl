using Kronecker
using LinearAlgebra

N = 3
# Note: assuming symmetric, Nx = Ny
# If using Nx and Ny instead of N, it is largely for readability / clarity
Nx, Ny = N, N
dx, dy = 1/N, 1/N
# if dx = dy, can abstract to h
h = dx
rho = 1
dxi, dyi, hi, dti, rhoi = 1/dx, 1/dy, 1/h, 1/dt, 1/rho

# velocity BC for lid driven flow problem
u_top = 1
u_bottom, v_left, v_right = 0, 0, 0

# timestep chosen with CFL condition uΔt/Δx < 1
dt = 0.1*h

# NOT entirely original implementation here, testing guide
# what about boundary conditions, sentinel?
# the plus one might be handling? unclear

# convenienet indeces when iterating
# BC prevent starting count from 1
imin, jmin = 2, 2
imax, jmax = Nx + 1, Ny + 1

# interleave u and v update to mitigate cache misses for large matrices?
# indexing is flipped here!
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
BC = zeros(N, N)
BC[1, 1] = -1
BC[N, N] = -1
P += (Id ⊗ BC)
P[1:N, 1:N] += -1*Id
P[N*(N-1)+1:N^2, N*(N-1)+1:N^2] += -1*Id
P[1, :] .= 0
P[1, 1] = 1

# RHS of Pressure Poisson Equation
# do not reinitialize R every single time
# indeces are flipped!
# to fill out in column major order for both keep u transposed?
#   can check if makes differences
R = zeros(N^2)
u = 1
v = [1; 2]
for j in jmin:jmax
    for i in imin:imax
        n += 1
        R[n] = -rho*dti*((u[i+1, j] - u[i, j])*dxi + (v[i, j+1] - v[i, j])*dyi)
    end
end
    
# example using a direct solve
p = P \ R
# use a view+reshape to use pressure in matrix form
# solver will output vector
p_matrix = @view p
p_matrix = reshape(p_matrix, Ny, Nx)

# velocity correction with updated pressure
# flipped!
# tranpose or fuse loops somehow?
# also iterating from jmin to jmax may be super weird?
# if jmin is on bottom then traversing up the rows?
# can benchmark against linear solve and see how consequential
# use SIMD, @SIMD?
for j in jmin:jmax
    for i=imin+1:imax
        u[i, j] -= dt*rhoi*dxi* (p[i, j] - p[i-1, j])
    end
end
for j=jmin+1: jmax
    for i=imin : imax
        v[i, j] -= dt*rhoi*dyi* (p[i, j] - p[i, j-1])
    end
end
    
# BC
# BC for pressure already handled in the modified 2D Laplacian matrix
# BC for velocity using fictitious velocities
# perhaps can fuse these into previous loops so don't have huge cache misses
# just pad the initial matrices for v and u with 1 or 2 layers? (depending on whether us jmax+1 and jmin-1, etc)
# flipped
u[:, jmin-1] = u[:, jmin] - 2*(u[:, jmin] - u_bottom)
u[:, jmax+1] = u[:, jmax] - 2*(u[:, jmax] - u_top)
v[imin-1, :] = v[imin, :] - 2*(v[imin, :] - v_left)
v[imax+1, :] = v[imax, :] - 2*(v[imax, :] - v_right)

