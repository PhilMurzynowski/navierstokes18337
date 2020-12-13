
using Kronecker
using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap

"""
Maintain BC for velocity using ficitious velocities.
"""
function updateVelBC(u, v, vel_BC, opts)
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
    u[jmin-1, :] = u_jmin .- 2*(u_jmin .- u_top)
    u[jmax+1, :] = u_jmax .- 2*(u_jmax .- u_bottom)
    v[:, imin-1] = v_imin .- 2*(v_imin .- v_left)
    v[:, imax+1] = v_imax .- 2*(v_imax .- v_right)
    return
end

# NOTE! the signs may be incorrect!!
# double check
function predictVelSIMPLE(u, v, p_matrix, u_new, v_new, opts)
    Δt = opts["dt"]
    Δx, Δy = opts["dx"], opts["dy"]
    Re_inv = opts["Rei"]

    # work out bounds
    for i in imin+1:imax
        for j in jmin:jmax
            # linearly interpolate
            #v_interpolated = 1/4*(v[j, i] + v[j, i+1] + v[j-1, i] + v[j-1, i+1])
            v_interpolated = 0.25 * ( v[j, i-1] + v[j+1, i-1] + v[j, i] + v[j+1, i])
            # group together in stages for readability
            a3 = (u[j, i+1] - 2*u[j, i] + u[j, i-1]) / Δx^2
            a4 = (u[j+1, i] - 2*u[j, i] + u[j-1, i]) / Δy^2
            # calculations using interpolated velocites
            #a1 = -(u[j, i+1]^2 - u[j, i-1]^2)/(2*Δx) - v_interpolated*(u[j+1, i] - u[j-1, i])/(2*Δy)
            # try not squaring u, may have problems with BC?
            a1 = -(u[j, i])*(u[j, i+1] - u[j, i-1])/(2*Δx) - v_interpolated*(u[j+1, i] - u[j-1, i])/(2*Δy)
            # may be off by a sign here, a1, check
            # checked, looks fine
            A = a1 + Re_inv*(a3 + a4)
            # u^{n+1}
            #                                   edited indices here
            # p_matrix subtracts, j-1, i-1, as of different size, N x N
            u_new[j, i] = u[j, i] + Δt*(A - (p_matrix[j-1, i-1] - p_matrix[j-1, i-2])/Δx)
        end
    end
    for i in imin:imax
        for j in jmin+1:jmax
            # linearly interpolate
            #u_interpolated = 1/4*(u[j, i] + u[j, i-1] + u[j+1, i] + u[j+1, i-1])
            u_interpolated = 1/4*(u[j, i+1] + u[j, i] + u[j-1, i+1] + u[j-1, i])
            # group together in stages for readability
            b3 = (v[j+1, i] - 2*v[j, i] + v[j-1, i]) / Δy^2
            b4 = (v[j, i+1] - 2*v[j, i] + v[j, i-1]) / Δx^2
            # calculations using interpolated velocites
            # try not squaring v? BC problems with fictitious velocity?
            #b1 = -(v[j+1, i]^2 - v[j-1, i]^2)/(2*Δy) - u_interpolated*(v[j, i+1] - v[j, i-1])/(2*Δx)
            b1 = -v[j, i]*(v[j+1, i] - v[j-1, i])/(2*Δy) - u_interpolated*(v[j, i+1] - v[j, i-1])/(2*Δx)
            # may be off by a sign here, b1, check
            # checked, I think postive b1 is actually correct
            B = b1 + Re_inv*(b3 + b4)
            # v^{n+1}
            #                                   edited indices here
            # p_matrix subtracts, j-1, i-1, as of different size, N x N
            v_new[j, i] = v[j, i] + Δt*(B - (p_matrix[j-1, i-1] - p_matrix[j-2, i-1])/Δy)
        end
    end
    return
end

#=
saving code snippet
    # linearly interpolate
    u_interpolated = 1/4*(u[j, i] + u[j, i-1] + u[j+1, i] + u[j+1, i-1])
    v_interpolated = 1/4*(v[j, i] + v[j, i+1] + v[j-1, i] + v[j-1, i+1])
    # group together in stages for readability
    a3 = (u[j, i+1] - 2*u[j, i] + u[j, i-1]) / Δx^2
    a4 = (u[j+1, i] - 2*u[j, i] + u[j-1, i]) / Δy^2
    b3 = (v[j+1, i] - 2*v[j, i] + v[j-1, i]) / Δy^2
    b4 = (v[j, i+1] - 2*v[j, i] + v[j, i-1]) / Δx^2
    # calculations using interpolated velocites
    a1 = -(u[j, i+1]^2 - u[j, i-1]^2)/(2*Δx) - v_interpolated*(u[j+1, i] - u[j-1, i])/(2*Δy)
    b1 = -(v[j+1, i]^2 - v[j-1, i]^2)/(2*Δy) - u_interpolated*(v[j, i+1] - v[j, i-1])/(2*Δx)
    # believe paper wrongly negated a1 here, check
    A = a1 + Re_inv*(a3 + a4)
    B = b1 + Re_inv*(b3 + b4)
    # u^{n+1}
    u_new[j, i] = u[j, i] + Δt*(A - (p[j, i+1] - p[j, i])/Δx)
    v_new[j, i] = v[j, i] + Δt*(B - (p[j+1, i] - p[j, i])/Δy)
=#

# RHS of Pressure Poisson Equation
# to fill out in column major order for both keep u transposed?
#   can check if makes differences
function poissonRSIMPLE(u, v, R, opts)
    jmin = opts["jmin"]
    jmax = opts["jmax"]
    imin = opts["imin"]
    imax = opts["imax"]
    rho = opts["rho"]
    dti = opts["dti"]
    dxi = opts["dxi"]
    dyi = opts["dyi"]

    k = 1
    # j must be inner loop if using column major order
    for i in imin:imax
        for j in jmin:jmax
            #R[k] = -rho*dti*((u[j, i+1] - u[j, i])*dxi + (v[j+1, i] - v[j, i])*dyi)
            # different papers disagree on whether negated or not
            R[k] = -dti*((u[j, i+1] - u[j, i])*dxi + (v[j+1, i] - v[j, i])*dyi)
            #R[k] = dti*((u[j, i+1] - u[j, i])*dxi + (v[j+1, i] - v[j, i])*dyi)
            # switched indices
            # R[k] = -dti*((u[j, i] - u[j, i-1])*dxi + (v[j, i] - v[j-1, i])*dyi)
            k += 1
        end
    end
end

function pressureSolveSIMPLE(P, R, opts)
    # direct solve currently
    p = P \ R
    return p 
end

# currently updating u_corrector, v_corrector
# perhaps can modify this later add to u, v directly?
#   would that have complications?
# u and v have offsets [j+1, i+1] because p_corrector_mtx is only Ny*Nx
function velCorSIMPLE(u, v, p_corrector_mtx, opts)
    Δt = opts["dt"]
    Δx, Δy = opts["dx"], opts["dy"]
    Ny = opts["Ny"]
    Nx = opts["Nx"]
    # work out exact bounds
    # fuse loops, tranpose?
    for i in 2:Nx
        for j in 1:Ny
            # note this must be minus, as u' is -1*RHS of eq below
            u[j+1, i+1] -= 1 / (Δt*Δx) * (p_corrector_mtx[j, i] - p_corrector_mtx[j, i-1])
        end
    end
    for i in 1:Nx
        for j in 2:Ny
            v[j, i] -= 1 / (Δt*Δy) * (p_corrector_mtx[j, i] - p_corrector_mtx[j-1, i])
        end
    end
    return
end

function plot_uvp(u, v, p, opts)
    h = opts["h"]
    N = opts["N"]
   
    # have to do some transposing and reversing here for formatting
    # reverse ys because jmin is near top, while origin is bottom left
    #display(u)
    #display(v)
    # using quiver
    # not sure if can use streamplot
    xs = 0.0:h:h*(N+1)
    ys = reverse(0.0:h:h*(N+1))
    quiv = quiver(xs, ys, u', v', arrowsize = 0.1)

    # pressure
    # testing
    #p .= 0
    #p[1] = 1
    #p[2] = 2
    #p[3] = 3
    #p[65] = 4
    # transpose when reshaping here??
    # might need pressure reshaped into row major because formed laplacian that way
    pressure = reshape(p, N, N)'
    display(pressure)
    xs = 0.0:h:h*(N-1)
    ys = 0.0:h:h*(N-1)
    #ys = reverse(0.0:h:h*(N-1))
    hm = heatmap(xs, ys, pressure)
    cl = colorlegend(hm[end], raw = true, camera = campixel!)

    parent = Scene(resolution= (1000, 500))
    full_scene = vbox(vbox(hm, cl), quiv, parent=parent)
    display(full_scene)
    return
end

function BCLaplacianSIMPLE(opts)
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
    P[1, 1] = 1*hi^2
    #P[N^2, :] .= 0
    #P[N^2, 1] = 1*hi^2
    return P
end

function runSIMPLE(opts)
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
    p_corrector = zeros(Ny*Nx)
    # R is RHS of Pressure Poisson Equation
    R = zeros(Ny*Nx)
    # Laplacian with BC
    P = BCLaplacianSIMPLE(opts)

    # Velocity BC
    # velocity BC for lid driven flow problem
    vel_BC = Dict("u_top"=>0.5,
                  "u_bottom"=>0.0,
                  "v_left"=>0.0,
                  "v_right"=>0.0)

    t0, T = opts["t0"], opts["T"]
    # time loop
    for t in t0:dt:T
        updateVelBC(u, v, vel_BC, opts)
        p_matrix = @view p[:, :]
        p_matrix = reshape(p_matrix, Ny, Nx)
        predictVelSIMPLE(u, v, p_matrix, u_new, v_new, opts)
        u, u_new = u_new, u # swap, memory reuse
        v, v_new = v_new, v
        #println("u")
        #display(u)
        #println("v")
        #display(v)
        #poissonRSIMPLE(u, v, R, opts)
        poissonRSIMPLE(u, v, R, opts)
        #println("R")
        #display(R)
        # debugging
        p_corrector = pressureSolveSIMPLE(P, R, opts)
        # view and reshape here for readability
        # is this reshaping in an unexpected manner?
        # transpose or not?
        p_corrector_mtx = reshape(p_corrector, Ny, Nx)
        #println("p_corrector_mtx")
        #display(p_corrector_mtx)
        #println("before")
        #display(u)
        #display(v)
        velCorSIMPLE(u, v, p_corrector_mtx, opts)
        #println("after")
        #display(u)
        #display(v)
        #println("corrected u")
        #display(u)
        #println("corrected v")
        #display(v)
        p += p_corrector
        #plot_uvp(u, v, p, opts)
        #return
    end

    # plot
    # remember outermost boundaries are fictitious!
    display(u)
    display(v)
    #display(p_corrector)
    plot_uvp(u, v, p, opts)
    return
end



# highlevel parameters
# number of cells in one axis
# Note: currenlty assuming symmetric, Nx = Ny
# If using Nx and Ny instead of N, it is largely for readability / clarity
N = 8
Nx, Ny = N, N
h = 1
#h = 1/100
# if dx = dy, can abstract to h
dx, dy = h, h
# fluid parameters
rho = 100   # may not be used depending on impl
mu = .1     # may not be used depending on impl
Re = 100    # may not be used depending on impl
# timestep chosen with CFL condition uΔt/Δx < 1
dt = 1
# inverses for convenience
dxi, dyi, hi, dti, rhoi, Rei = 1/dx, 1/dy, 1/h, 1/dt, 1/rho, 1/Re
# convenienet indeces when iterating
# BC prevent starting count 
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
            "rhoi"=>rhoi,
            "mu"=>mu,
            "Re"=>Re,
            "Rei"=>Rei,
            "imin"=>imin,
            "jmin"=>jmin,
            "imax"=>imax,
            "jmax"=>jmax,
            "t0"=>0.0,
            "T"=>20.0 # working for one timestep?
            )


runSIMPLE(opts)