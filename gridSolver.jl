using Kronecker
using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap

include("CG.jl")

function set_vel_BC!(u, v, opts_BC)
    u_top = opts_BC["u_top"]
    u_bottom = opts_BC["u_bottom"]
    v_left = opts_BC["v_left"]
    v_right = opts_BC["v_right"]

    # creating copies for debuggin
    #u_new = deepcopy(u)
    #v_new = deepcopy(v)


    # rest should be 0 if at BC
    # may be able to remove
    u[:, 1] .= 0
    u[:, end] .= 0
    v[1, :] .= 0
    v[end, :] .= 0

    #println("u_top")
    #println(u_top)
    u[1, :] .= u_bottom
    u[end, :] .= u_top
    v[:, 1] .= v_left
    v[:, end] .= v_right

    return
end

#@inline function velocity_update!(u1, v1, u2, v2, p, opts)
#function velocity_update(u1, v1, u2, v2, p, opts)
function velocity_update(u1, v1, p, opts)
    # inlining because many of these functions will be reusing views
    # use two arrays with swapping for memory reuse, u2, v2 can initially be garbage
    Δt = opts["dt"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    μ = opts["mu"]
    ρ = opts["rho"]

    u2 = deepcopy(u1)
    v2 = deepcopy(v1)

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

    #println("end of velocity update")
    #println("u2")
    #display(u2)
    #println("v2_c")
    #display(v2_c)
    # swapping inside a function won't do anything
    #v1, v2 = v2, v1
    #u1, u2 = u2, u1
    return u2, v2
end

"""
function to assemble term describing pressure dependence on velocity
    output: v_dependence
        passed in for memory reuse
"""
#function v2p(v_dependence, u, v, opts)
function v2p(u, v, opts)
    Δt = opts["dt"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    N = opts["N"]

    v_dependence = zeros(N, N)

    #println("inside v2p")
    #println("u")
    #display(u)
    #println("v")
    #display(v)

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

    #depend_inner = @view v_dependence[2:end-1, 2:end-1]
    # perhaps can do common subexpression elemination
    # lengthy and in one line to try to get compiler to optimize subexpression
    #depend_inner = (1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))



    v_dependence = (1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))
                    - 1/(2*Δx)^2 * (u_r - u_l).^2
                    - 1/(2*Δy)^2 * (v_b - v_t).^2 
                    - 1/(2*Δy)^2 * (v_b - v_t).^2
                    - (1/(Δy)*(u_b - u_t)) .* (1/(Δx)*(v_r - v_l)) * 1/2
                    )
    #depend_inner = 1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))
    #depend_inner -= 1/(2*Δx)^2 * (u_r - u_l).^2
    #depend_inner -= 1/(2*Δy)^2 * (v_b - v_t).^2
    #depend_inner -= (1/(2Δy)*(u_b - u_t)) .* (1/(2*Δx)*(v_r - v_l))


    return v_dependence
end


# test of v2p
#u = ones(N, N)*3
#v = ones(N, N)*2
#v_dependence = zeros(N, N)
#v2p(u, v, v_dependence, opts)
#v_dependence

#function grid_poisson_solve(p1, p2, velocity_component, opts)
    # use two arrays with swapping for memory reuse, p2 can initially be garbage
    # currently creating tmp manually for debugging
function grid_poisson_solve(p, velocity_component, opts)
    # currently using broadcast operations on matrices
    # iterating until convergence

    # maybe copy is enough?
    p1 = deepcopy(p)
    
    Δx = opts["dx"]
    Δy = opts["dy"]
    ρ = opts["rho"]
    iter = opts["poisson_iter"]

    # instead of using a while loop and comparing size of update in matrix,
    # simply run a number of iterations and ouput difference from last update
    # for analysis
    for i in 1:iter
        # maybe copy is enough?
        p1 = deepcopy(p)

        p1_r = @view p1[2:end-1, 3:end]
        p1_l = @view p1[2:end-1, 1:end-2]
        p1_b = @view p1[3:end, 2:end-1]
        p1_t = @view p1[1:end-2, 2:end-1]

        # can compute prefix to its own constant if not done by compiler
        p[2:end-1, 2:end-1] = (1.0/(2*(Δx^2 + Δy^2)) * ((Δx^2*(p1_t + p1_b)) + (Δy^2*(p1_r - p1_l)))
                            - ρ*Δx^2*Δy^2/(2*(Δx^2 + Δy^2)).*velocity_component
                            )

        # update BC
        p[:, 1] = p[:, 2]
        p[1, :] = p[2, :]
        p[:, end] = p[:, end-1]
        p[end, :] .= 0 # do this last to avoid overwriting

        # swap to reuse memory
        #p1, p2 = p2, p1
    end

    #return maximum(abs.(p2 - p1)), p1
    return maximum(abs.(p - p1)), p
    # desired answer will be in p1
end

# test of grid_poisson_solve
#u = ones(N, N)*3
#v = ones(N, N)*2
#v_dependence = zeros(N, N)
#v2p(u, v, v_dependence, opts)
#pressure1 = zeros(N+2, N+2)
#pressure2 = zeros(N+2, N+2)
#pressure_eps = grid_poisson_solve(pressure1, pressure2, v_dependence, opts)



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

# full test

function run_simulation(opts, BC_opts)
    
    sim_iter = opts["simulation_iter"]
    N = opts["N"]

    u1_mtx, v1_mtx = zeros(N+2, N+2), zeros(N+2, N+2)
    #u2_mtx, v2_mtx = zeros(N+2, N+2), zeros(N+2, N+2)
    #pressure1, pressure2 = zeros(N+2, N+2), zeros(N+2, N+2)
    pressure = zeros(N+2, N+2)
    v_for_poisson = zeros(N, N)

    set_vel_BC!(u1_mtx, v1_mtx, BC_opts)

    for s in 1:sim_iter
        #v_for_poisson = v2p(v_for_poisson, u1_mtx, v1_mtx, opts)
        v_for_poisson = v2p(u1_mtx, v1_mtx, opts)
        #println("in sim loop")
        #println("v_dependence")
        #display(v_dependence)
        #pressure_eps, pressure1 = grid_poisson_solve(pressure1, pressure2, v_for_poisson, opts)
        pressure_eps, pressure = grid_poisson_solve(pressure, v_for_poisson, opts)
        #println(pressure_eps) # debuggin
        #u1_mtx, v1_mtx = velocity_update(u1_mtx, v1_mtx, u2_mtx, v2_mtx, pressure1, opts)
        #u1_mtx, v1_mtx = velocity_update(u1_mtx, v1_mtx, u2_mtx, v2_mtx, pressure, opts)
        u1_mtx, v1_mtx = velocity_update(u1_mtx, v1_mtx, pressure, opts)
        #set_vel_BC!(u1_mtx, v1_mtx, BC_opts)
        #println("u1")
        #display(u1)
        #println("v1")
        #display(v1)
        #if s == 4
        #    #plot_uvp(u1_mtx, v1_mtx, pressure1, opts)
        #    plot_uvp(u1_mtx, v1_mtx, pressure, opts)
        #    return
        #end

        #set_vel_BC(u1, v1, BC_opts)
        #set_vel_BC(u2, v2, BC_opts)
    end
    
    #plot_uvp(u1_mtx, v1_mtx, pressure1, opts)
    plot_uvp(u1_mtx, v1_mtx, pressure, opts)
    #display(pressure1)
    #display(u1)
    #display(v1)
    #println("dual")
    #display(pressure2)
    #display(u2)
    #display(v2)

end

# can declare as constants

N = 32
h = 3/N # chose this spacing largely for visuals
dt = 0.001
rho = 10.0
mu = 0.15
opts = Dict("poisson_iter"=>50,
            "simulation_iter"=>500,
            "N"=>N,
            "Nx"=>N,
            "Ny"=>N,
            "dx"=>h,
            "dy"=>h,
            #"dxi"=>dxi,
            #"dyi"=>dyi,
            "h"=>h,
            #"hi"=>hi,
            "dt"=>dt,
            #"dti"=>dti,
            "rho"=>rho,
            "mu"=>mu,
            #"rhoi"=>rhoi,
            #"mu"=>mu,
            #"Re"=>Re,
            #"Rei"=>Rei,
            #"imin"=>imin,
            #"jmin"=>jmin,
            #"imax"=>imax,
            #"jmax"=>jmax,
            #"t0"=>0.0,
            #"T"=>0.3 # working for one timestep?
            )
BC_opts = Dict("u_top"=>1.0,
               "u_bottom"=>0.0,
               "v_left"=>0.0,
               "v_right"=>0.0)


run_simulation(opts, BC_opts)