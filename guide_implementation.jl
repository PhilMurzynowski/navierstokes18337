using Kronecker
using LinearAlgebra

function predict_vel(u, v, u_new, v_new, opts)
    jmin = opts["jmin"]
    jmax = opts["jmax"]
    imin = opts["imin"]
    imax = opts["imax"]
    dt = opts["dt"]
    dxi = opts["dxi"]
    dyi = opts["dyi"]
    mu = opts["mu"]
    # interleave u and v update to mitigate cache misses for large matrices?
    for j in jmin+1:jmax
        for i in imin:imax
            u_interpolated = 0.25 * ( u[j-1, i] + u[j, i] + u[j-1, i+1] + u[j, i+1])
            v_update = (v[j, i-1] - 2*v[j, i] + v[j, i+1]) * dxi^2 * mu
            v_update += (v[j-1, i] - 2*v[j, i] + v[j+1, i]) * dyi^2 * mu
            v_update -= v[j, i] * (v[j+1, i] - v[j-1, i]) * 0.5 * dyi
            v_update -= u_interpolated * (v[j, i+1] - v[j, i-1]) * 0.5 * dxi
            v_update *= dt
            v_new[j, i] = v[j, i] + v_update
        end
    end
    for j in jmin:jmax
        for i in imin+1:imax
            v_interpolated = 0.25 * ( v[j, i-1] + v[j+1, i-1] + v[j, i] + v[j+1, i])
            u_update = (u[j, i-1] - 2*u[j, i] + u[j, i+1]) * dxi^2 * mu
            u_update += (u[j-1, i] - 2*u[j, i] + u[j+1, i]) * dyi^2 * mu
            u_update -= u[j, i] * (u[j, i+1] - u[j, i-1]) * 0.5 * dxi
            u_update -= v_interpolated * (u[j+1, i] - u[j-1, i]) * 0.5 * dyi
            u_update *= dt
            u_new[j, i] = u[j, i] + u_update
        end
    end
end

# 2D Laplacian with BC
# can describe the general structure with Kronecker products
# with some small modifications due to Von Neumman boundary conditions
# and also first cell is kept constant as pressure is defined up to a constant
#   other cells in reference to cell that is kept constant
# L is tridiagonal Laplacian Matrix
function init_BCLaplacian(opts)
    N = opts["N"]
    hi = opts["hi"]

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
    return P
end

# RHS of Pressure Poisson Equation
# to fill out in column major order for both keep u transposed?
#   can check if makes differences
function update_poisson_RHS(u, v, R, opts)
    jmin = opts["jmin"]
    jmax = opts["jmax"]
    imin = opts["imin"]
    imax = opts["imax"]
    rho = opts["rho"]
    dti = opts["dti"]
    dxi = opts["dxi"]
    dyi = opts["dyi"]

    k = 1
    for j in jmin:jmax
        for i in imin:imax
            R[k] = -rho*dti*((u[j, i+1] - u[j, i])*dxi + (v[j+1, i] - v[j, i])*dyi)
            k += 1
        end
    end
end

function pressure_solve(P, R, p, opts)
    # direct solve currently
    p = P \ R
end

# velocity correction with updated pressure
# tranpose or fuse loops somehow?
# also iterating from jmin to jmax may be super weird?
# if jmin is on bottom then traversing up the rows?
# can benchmark against linear solve and see how consequential
# use SIMD, @SIMD?
# u and v have offsets [j+1, i+1] because p_matrix is only Ny*Nx
function correct_vel(u, v, p_matrix, opts)
    Ny = opts["Ny"]
    Nx = opts["Nx"]
    dt = opts["dt"]
    rhoi = opts["rhoi"]
    dxi = opts["dxi"]
    dyi = opts["dyi"]

    for j in 1:Ny
        for i in 2:Nx
            u[j+1, i+1] -= dt*rhoi*dxi* (p_matrix[j, i] - p_matrix[j, i-1])
        end
    end
    for j in 2:Ny
        for i in 1:Nx
            v[j+1, i+1] -= dt*rhoi*dyi* (p_matrix[j, i] - p_matrix[j-1, i])
        end
    end
end
    
# BC
# BC for pressure already handled in the modified 2D Laplacian matrix
# BC for velocity using fictitious velocities
# perhaps can fuse these into previous loops so don't have huge cache misses
# just pad the initial matrices for v and u with 1 or 2 layers? (depending on whether us jmax+1 and jmin-1, etc)
# fictitious velocity approach seems superior as get exactly desired velocity on the boundary of the cell
# for now at least using views to not allocate more memory
#u[jmin-1, :] = u[jmin, :] - 2*(u[jmin, :] - u_bottom)
#u[jmax+1, :] = u[jmax, :] - 2*(u[jmax, :] - u_top)
#v[:, imin-1] = v[:, imin] - 2*(v[:, imin] - v_left)
#v[:, imax+1] = v[:, imax] - 2*(v[:, imax] - v_right)

function update_vel_BC(u, v, vel_BC, opts)
    jmin = opts["jmin"]
    jmax = opts["jmax"]
    imin = opts["imin"]
    imax = opts["imax"]

    u_top = vel_BC["u_top"]
    u_bottom = vel_BC["u_bottom"]
    v_left = vel_BC["v_left"]
    v_right = vel_BC["v_right"]

    u_jmin = @view u[jmin, :]
    u_jmax = @view u[jmax, :]
    v_imin = @view v[:, imin]
    v_imax = @view v[:, imax]
    u[jmin-1, :] = u_jmin .- 2*(u_jmin .- u_bottom)
    u[jmax+1, :] = u_jmax .- 2*(u_jmax .- u_top)
    v[:, imin-1] = v_imin .- 2*(v_imin .- v_left)
    v[:, imax+1] = v_imax .- 2*(v_imax .- v_right)
end



function runExample(opts)
    Ny = opts["Ny"]
    Nx = opts["Nx"]
    # can make matrices static size as volume unchanging
    # initialize all arrays with initial conditions
    #   pressure (p), u, v zero everywhere
    u, v = zeros(Ny+2, Nx+2), zeros(Ny+2, Nx+2)
    u_new, v_new = zeros(Ny+2, Nx+2), zeros(Ny+2, Nx+2)
    # p is the pressure in the region of interest
    # will be the solution to solving the poisson equation solve
    p = zeros(Ny*Nx)
    # R is RHS of Pressure Poisson Equation
    R = zeros(Ny*Nx)
    p_matrix = @view p[:, :]
    p_matrix = reshape(p_matrix, Ny, Nx)
    # Laplacian with BC
    P = init_BCLaplacian(opts)

    # Velocity BC
    # velocity BC for lid driven flow problem
    vel_BC = Dict("u_top"=>1,
                  "u_bottom"=>0,
                  "v_left"=>0,
                  "v_right"=>0)

    t0, T = opts["t0"], opts["T"]
    # time loop
    for t in t0:dt:T
        update_vel_BC(u, v, vel_BC, opts)
        predict_vel(u, v, u_new, v_new, opts)
        u, u_new = u_new, u # swap, memory reuse
        v, v_new = v_new, v
        #println("u")
        #display(u)
        #println("v")
        #display(v)
        #println("u_new")
        #display(u_new)
        #println("v_new")
        #display(v_new)
        update_poisson_RHS(u, v, R, opts)
        pressure_solve(P, R, p, opts)
        correct_vel(u, v, p_matrix, opts)
        # can option to plot here
    end
    # plot
    # they come out upside down! fix!
    # remember outermost boundaries are fictitious!
    display(u)
    display(v)
    display(p)
end

# NOT entirely original implementation here, testing guide

# initialize problem
# highlevel parameters
# number of cells in one axis
# Note: currenlty assuming symmetric, Nx = Ny
# If using Nx and Ny instead of N, it is largely for readability / clarity
N = 5
Nx, Ny = N, N
dx, dy = 1/N, 1/N
h = dx # if dx = dy, can abstract to h
# fluid parameters
rho = 1000
mu = 1
# timestep chosen with CFL condition uΔt/Δx < 1
dt = 0.1*h
# inverses for convenience
dxi, dyi, hi, dti, rhoi = 1/dx, 1/dy, 1/h, 1/dt, 1/rho
# convenienet indeces when iterating
# BC prevent starting count from 1
imin, jmin = 2, 2
imax, jmax = Nx + 1, Ny + 1

# group params
opts = Dict("N"=>N,
            "Nx"=>Nx,
            "Ny"=>Ny,
            "dx"=>dx,
            "dy"=>dy,
            "dxi"=>dxi,
            "dyi"=>dyi,
            "h"=>h,
            "hi"=>hi,
            "dt"=>dt,
            "dti"=>dti,
            "rho"=>rho,
            "mu"=>mu,
            "rhoi"=>rhoi,
            "imin"=>imin,
            "jmin"=>jmin,
            "imax"=>imax,
            "jmax"=>jmax,
            "t0"=>0.0,
            "T"=>0.2
            )


runExample(opts)