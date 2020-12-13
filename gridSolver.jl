
"""
function to assemble term describing pressure dependence on velocity
    output: v_dependence
        passed in for memory reuse
"""
function v2p(u, v, v_dependence, opts)
    Δt = opts["dt"]
    Δx = opts["dx"]
    Δy = opts["dy"]

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

    depend_inner = @view v_dependence[2:end-1, 2:end-1]
    # perhaps can do common subexpression elemination
    # lengthy and in one line to try to get compiler to optimize subexpression
    depend_inner = (1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))
                    - 1/(2*Δx)^2 * (u_r - u_l).^2
                    - 1/(2*Δy)^2 * (v_b - v_t).^2 
                    - 1/(2*Δy)^2 * (v_b - v_t).^2
                    - (1/(2Δy)*(u_b - u_t)) .* (1/(2*Δx)*(v_r - v_l))
                    )
    #depend_inner = 1/Δt * ((1/(2*Δx)*(u_r - u_l)) + 1/(2*Δy)*(v_b - v_t))
    #depend_inner -= 1/(2*Δx)^2 * (u_r - u_l).^2
    #depend_inner -= 1/(2*Δy)^2 * (v_b - v_t).^2
    #depend_inner -= (1/(2Δy)*(u_b - u_t)) .* (1/(2*Δx)*(v_r - v_l))
end

# can declare as constants
N = 7
h = 1.0
dt = 1.0
rho = 1.0
opts = Dict("N"=>N,
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
            #"rhoi"=>rhoi,
            #"mu"=>mu,
            #"Re"=>Re,
            #"Rei"=>Rei,
            #"imin"=>imin,
            #"jmin"=>jmin,
            #"imax"=>imax,
            #"jmax"=>jmax,
            "t0"=>0.0,
            "T"=>0.3 # working for one timestep?
            )

# test of v2p
#u = ones(N, N)*3
#v = ones(N, N)*2
#v_dependence = zeros(N, N)
#v2p(u, v, v_dependence, opts)
#v_dependence

function poisson_solve(p1, p2, opts)
    # currently using broadcast operations on matrices
    # iterating until convergence
    # use two arrays with swapping for memory reuse

    p1_r = @p1iew p1[2:end-1, 3:end]
    p1_l = @p1iew p1[2:end-1, 1:end-2]
    p1_b = @p1iew p1[3:end, 2:end-1]
    p1_t = @p1iew p1[1:end-2, 2:end-1]

    p2 = (Δx^2*(p1_t + p1_b))+ (Δy^2*(p1_r - p1_l))


    # return number of iterations for analysis
    return iter
end