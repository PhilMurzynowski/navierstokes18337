using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap
using Printf

# file name not capitalized to GridSolver.jl by accident

"""
Few operations to set BC
Inlined as only a few slicing operations
"""
@inline function set_vel_BC!(u, v, opts_BC)
    u_top = opts_BC["u_top"]
    u_bottom = opts_BC["u_bottom"]
    v_left = opts_BC["v_left"]
    v_right = opts_BC["v_right"]
    # rest e.g u_right should be 0 if at boundary as can't go through wall
    # don't actually need to update since init to 0 and only interior points updated
    #=
    u[:, 1] .= 0
    u[:, end] .= 0
    v[1, :] .= 0
    v[end, :] .= 0
    =#
    u[1, :] .= u_bottom
    u[end, :] .= u_top
    v[:, 1] .= v_left
    v[:, end] .= v_right

    return
end

"""
Function to update velocity
Inlining because many of these functions will be reusing views,
    perhaps could get picked up by compiler
Use two arrays with swapping for memory reuse, u2, v2 can initially be garbage
    in the inner points, the BC must be set properly for both
"""
@inline function velocity_update(u1, v1, u2, v2, p, opts)
    Δt = opts["dt"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    μ = opts["mu"]
    ρ = opts["rho"]

    u1_c = @view u1[2:end-1, 2:end-1]
    u1_r = @view u1[2:end-1, 3:end]
    u1_l = @view u1[2:end-1, 1:end-2]
    u1_b = @view u1[3:end, 2:end-1]
    u1_t = @view u1[1:end-2, 2:end-1]

    v1_c = @view v1[2:end-1, 2:end-1]
    v1_r = @view v1[2:end-1, 3:end]
    v1_l = @view v1[2:end-1, 1:end-2]
    v1_b = @view v1[3:end, 2:end-1]
    v1_t = @view v1[1:end-2, 2:end-1]

    p_r = @view p[2:end-1, 3:end]
    p_l = @view p[2:end-1, 1:end-2]
    p_b = @view p[3:end, 2:end-1]
    p_t = @view p[1:end-2, 2:end-1]

    u2[2:end-1, 2:end-1] = (u1_c +
                        Δt*(-1/Δx*u1_c.*(u1_c - u1_l)
                            - 1/Δy*v1_c.*(u1_c - u1_t)
                            - 1/(2*ρ*Δx)*(p_r - p_l)
                            + μ/Δx^2*(u1_r - 2*u1_c + u1_l)
                            + μ/Δy^2*(u1_b - 2*u1_c + u1_t)
                            ))

    v2[2:end-1, 2:end-1] = (v1_c +
                        Δt*(-1/Δy*v1_c.*(v1_c - v1_t)
                            - 1/Δx*u1_c.*(v1_c - v1_l)
                            - 1/(2*ρ*Δy)*(p_b - p_t)
                            + μ/Δx^2*(v1_r - 2*v1_c + v1_l)
                            + μ/Δy^2*(v1_b - 2*v1_c + v1_t)
                            ))

    return u2, v2
end

"""
function to assemble term describing pressure dependence on velocity
    output: v_dependence
        passed in for memory reuse
"""
function v2p(v_dependence, u, v, opts)
    Δt = opts["dt"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    N = opts["N"]

    # using views to attempt to minimize creating new arrays
    # using this indexing because have to be careful with BC, boundary conditions
    u_r = @view u[2:end-1, 3:end]
    u_l = @view u[2:end-1, 1:end-2]
    u_b = @view u[3:end, 2:end-1]
    u_t = @view u[1:end-2, 2:end-1]

    v_r = @view v[2:end-1, 3:end]
    v_l = @view v[2:end-1, 1:end-2]
    v_b = @view v[3:end, 2:end-1]
    v_t = @view v[1:end-2, 2:end-1]

    # lengthy and in one line to try to get compiler to optimize / eliminate subexpressions
    v_dependence = (1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))
                    - 1/(2*Δx)^2 * (u_r - u_l).^2
                    - 1/(2*Δy)^2 * (v_b - v_t).^2 
                    - 1/(2*Δy)^2 * (v_b - v_t).^2
                    - (1/(Δy)*(u_b - u_t)) .* (1/(Δx)*(v_r - v_l)) * 1/2
                    )

    return v_dependence
end

"""
Solve for pressure by continuially iterating over global array
    and using local updates until convergence.
Uses two arrays with swapping for memory reuse, p2 can initially be garbage
Currently using broadcast operations on matrices for vectorization
"""
function grid_poisson_solve(p1, p2, velocity_component, opts, max_iter=500)
    
    Δx = opts["dx"]
    Δy = opts["dy"]
    ρ = opts["rho"]
    ϵ = opts["ϵ"]

    iter = 0
    residual = Inf

    # in reality could check residual much more periodically or run
    # for a set number of iterations
    # as subtracting matrices and norm operation quite expensive to do every operation
    while residual > ϵ^2 && iter < max_iter
        iter += 1

        p1_r = @view p1[2:end-1, 3:end]
        p1_l = @view p1[2:end-1, 1:end-2]
        p1_b = @view p1[3:end, 2:end-1]
        p1_t = @view p1[1:end-2, 2:end-1]

        # can compute prefix to its own constant if not done by compiler
        p2[2:end-1, 2:end-1] = (1.0/(2*(Δx^2 + Δy^2)) * ((Δx^2*(p1_t + p1_b)) + (Δy^2*(p1_r - p1_l)))
                            - ρ*Δx^2*Δy^2/(2*(Δx^2 + Δy^2)).*velocity_component
                            )

        # update BC
        p2[:, 1] = p2[:, 2]
        p2[1, :] = p2[2, :]
        p2[:, end] = p2[:, end-1]
        p2[end, :] = p2[end-1, :]

        residual = norm(p2 - p1)
        p2, p1 = p1, p2
        #@printf "iter: %d, residual: %10f\n" iter residual
    end

    return p2, iter
end

function plot_uvp(u, v, p, opts)
    h = opts["h"]
    N = opts["N"]
   
    # have to do some transposing here for formatting
    # pretty sure not the cleanest way but got annoyed and left it for now
    xs = 0.0:h:h*(N+1)
    ys = 0.0:h:h*(N+1)
    quiv = quiver(xs, ys, u', v', arrowsize = 0.01)

    p_xs = 0.0:h:h*(N-1)
    p_ys = 0.0:h:h*(N-1)
    hm = AbstractPlotting.heatmap(p_xs, p_ys, p')
    cl = colorlegend(hm[end], raw = true, camera = campixel!)

    parent = Scene(resolution= (1000, 500))
    full_scene = vbox(vbox(hm, cl), quiv, parent=parent)
    display(full_scene)
    return
end


function run_simulation(opts, BC_opts)
    
    sim_iter = opts["simulation_iter"]
    N = opts["N"]

    # preallocate
    u1_mtx, v1_mtx = zeros(N+2, N+2), zeros(N+2, N+2)
    u2_mtx, v2_mtx = zeros(N+2, N+2), zeros(N+2, N+2)
    pressure1, pressure2 = zeros(N+2, N+2), zeros(N+2, N+2)
    v_for_poisson = zeros(N, N)

    # as will be swapping u1, u2 and v1, v2 make sure both have BC
    # BC will not be affected in future as only update inner points
    set_vel_BC!(u1_mtx, v1_mtx, BC_opts)
    set_vel_BC!(u2_mtx, v2_mtx, BC_opts)

    for s in 1:sim_iter
        v_for_poisson = v2p(v_for_poisson, u1_mtx, v1_mtx, opts)
        pressure1, num_iter = grid_poisson_solve(pressure1, pressure2, v_for_poisson, opts)
        #pressure1, num_iter = grid_poisson_solve(pressure1, v_for_poisson, opts)
        u1_mtx, v1_mtx = velocity_update(u1_mtx, v1_mtx, u2_mtx, v2_mtx, pressure1, opts)
    end
    
    return u1_mtx, v1_mtx, pressure1

end


N = 32
h = 3/N # chose this spacing largely for visuals
dt = 0.001
# have to make rho 1 as didn't get around to including rho in later iterations
# rho = 10.0
rho = 1.0
mu = 0.15
opts = Dict("simulation_iter"=>100,                
            "ϵ"=>1e-4,          
            "N"=>N,
            "Nx"=>N,
            "Ny"=>N,
            "dx"=>h,
            "dy"=>h,
            "h"=>h,
            "dt"=>dt,
            "rho"=>rho,
            "mu"=>mu,
            )
BC_opts = Dict("u_top"=>1.0,
               "u_bottom"=>0.0,
               "v_left"=>0.0,
               "v_right"=>0.0)


u, v, p = run_simulation(opts, BC_opts)
plot_uvp(u, v, p, opts)