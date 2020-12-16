using LinearAlgebra, Kronecker

include("CG.jl")

# general functions

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
    #BC[1, 1] = -hi^2
    #BC[N, N] = -hi^2
    #P += (Id ⊗ BC)
    #P[1:N, 1:N] += -hi^2*Id
    #P[N*(N-1)+1:N^2, N*(N-1)+1:N^2] += -hi^2*Id
    #P[1, :] .= 0
    #P[1, 1] = 1*hi^2
    #P[N^2, :] .= 0
    #P[N^2, 1] = 1*hi^2
    #return P
    return -P
end

# compute BC

# could this be fused to updateVoriticity for better cache use, repeating this question everywhere ha
function updateVorticityWallBC!(ω, ψ, opts, opts_BC)
    h = opts["h"]
    N = opts["N"]

    u_top = opts_BC["u_top"]
    u_bottom = opts_BC["u_bottom"]
    v_left = opts_BC["v_left"]
    v_right = opts_BC["v_right"]

    # boundary conditions are constant by problem specification

    #=
    c2 = @view ψ[:, 3]
    c3 = @view ψ[:, 2]
    cn = @view ψ[:, N]
    cn1 = @view ψ[:, N+1]

    r2 = @view ψ[2, :]
    r3 = @view ψ[3, :]
    rn = @view ψ[N, :]
    rn1 = @view ψ[N+1, :]

    # are views more efficient than iterating
    # or am I creating unnecessary temporary arrays
    # to specialize for lid driven cavity could also flip axes so
    #   that when updating top velocity only access a column
    ω[N+2, :] = 2/h^2 .* (rn - rn1) .- 2*u_top/h
    ω[1, :] = 2/h^2 .* (r3 - r2) .+ 2*u_bottom/h

    ω[:, 1] = 2/h^2 .* (c3 - c2) .- 2*v_left/h
    ω[:, N+2] = 2/h^2 .* (cn - cn1) .+ 2*v_right/h
    =#

    # iterating very likely better than slicing here
    # only have one term in BC as expecting other to be 0
    # bottom top will have lots of cache misses
    for i in 2:N-1
        # top
        # bottom
        ω[N, i] = -2/h^2*ψ[N-1, i] - 2/h*u_top
        ω[1, i] = -2/h^2*ψ[2, i] + 2/h*u_bottom
    end
    # left right
    for j in 2:N-1
        # left
        # right
        ω[j, 1] = -2/h^2*ψ[j, 2] + 2/h*v_left
        ω[j, N] = -2/h^2*ψ[j, N-1] - 2/h*v_right
    end
    return
end

# simulate vorticity forward a timestep
# starting with simple FTCS (forward time centered space) first order method
#   2nd order upwind for convection
#   check if worth over 1st order, not much harder to implement
# can switch to Dormand-Prince
# should I be storing u v in arrays
# does second order differencing have more potential for cache misses? consequence
function updateVorticity(ω, ω_next, ψ, opts)
    Δt = opts["dt"]
    h = opts["h"]   # just use h since know using dx and dy equal
    Re = opts["Re"]

    for i in 2:N-1
        #i_boundary = (i == 2 || i == N+1) ? true : false
        for j in 2:N-1
            # higher indices currently correspond to top
            # u = dψ/dy
            u_ji = 1/(2*h)*(ψ[j+1, i] - ψ[j-1, i])
            # v = -dψ/dx
            v_ji = -1/(2*h)*(ψ[j, i+1] - ψ[j, i-1])

            ω_next[j, i] = ω[j, i] + 
                        Δt*( -1/(2*h)*u_ji*(ω[j, i+1] - ω[j, i-1])
                            -1/(2*h)*v_ji*(ω[j+1, i] - ω[j-1, i])
                            + 1/(Re*h^2)*(ω[j, i+1] - ω[j, i-1] - 4*ω[j, i] + ω[j+1, i] + ω[j-1, i]))

            # if on the boundary don't use second order upwind
            # could use sentinels instead of control flow but more memory

            #= debugging, don't use 2nd order upwind at all for now

            if i_boundary || j == 2 || j == N + 1
                continue
            else
                ωx = (u_ji >= 0) ? 
                        (1/(3*h)*(-ω[j, i+2] + 3*ω[j, i+1] - 3*ω[j, i] + ω[j, i-1])) :
                        (1/(3*h)*(-ω[j, i+1] + 3*ω[j, i] - 3*ω[j, i-1] + ω[j, i-2]))
                ωy = (v_ji >= 0) ?
                        (1/(3*h)*(-ω[j+2, i] + 3*ω[j+1, i] - 3*ω[j, i] + ω[j-1, i])) :
                        (1/(3*h)*(-ω[j+1, i] + 3*ω[j, i] - 3*ω[j-1, i] + ω[j-2, i]))

                ω_next[j, i] += Δt*(
                             - 1/2*u_ji*ωx
                             - 1/2*v_ji*ωy
                             )
            end

            =#
        end
    end
    return ω_next
end


# solve poisson with CG

# is computing velocity necessary?
# believe can eliminate storing velocity


function runVorticityStream(opts, opts_BC)
    # init vorticity, streamfunction
    N = opts["N"]
    ω = zeros(N, N)
    ωtmp = zeros(N, N)
    ψ = zeros(N, N)
    # thought about it, and generating mtx of size other than NxN doesn't make sense
    # can remove param later
    P = genPoissonStreamMtx(N, opts)
    #println("P")
    #display(P)
    #return
    ϵ = 1e-10

    max_iter = 20

    for i in 1:max_iter
        println(i)
        updateVorticityWallBC!(ω, ψ, opts, opts_BC)
        println("BC")
        ω  = updateVorticity(ω, ωtmp, ψ, opts)
        println("update")
        #display(ω)
        #display(ψ)
        #return
        #ω_inner = @view ω[2:end-1, 2:end-1]
        #ω_vec = reshape(ω_inner, N*N, 1)
        ## sign
        #poisson_RHS = -ω_vec
        poisson_RHS = -1*reshape(ω, length(ω), 1)
        #poisson_RHS = ω_vec
        #ψ_inner = @view ψ[2:end-1, 2:end-1]
        #ψ_vec = vec(ψ_inner)

        #=
        ψ_vec = reshape(ψ, length(ψ), 1)
        println("enter CG")
        ψ_vec, num_iter = CG_std(P, poisson_RHS, ψ_vec, ϵ, N^2)
        #ψ[2:end-1, 2:end-1] = reshape(ψ_vec, N, N)
        =#
        ψ_vec = P \ poisson_RHS
        ψ = reshape(ψ_vec, N, N)

        println("leave CG")
    end

    display(ω)
    display(ψ)

end

N = 10
h = 1/N
dt = 0.0001
BC_opts = Dict("u_top"=>1.0,
               "u_bottom"=>0.0,
               "v_left"=>0.0,
               "v_right"=>0.0)
opts = Dict("N"=>N,
            "Nx"=>N,
            "Ny"=>N,
            "dx"=>h,
            "dy"=>h,
            "h"=>h,
            "dt"=>dt,
            "Re"=>150.0
            )
runVorticityStream(opts, BC_opts)