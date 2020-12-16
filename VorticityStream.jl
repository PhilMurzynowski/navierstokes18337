


# init

# general functions

function genPoissonMtx(opts)

    N = opts["N"]
    hi = 1/opts["h"]
    println("updated gen mtx")

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
    BC[1, 1] = -hi^2
    BC[N, N] = -hi^2
    P += (Id ⊗ BC)
    P[1:N, 1:N] += -hi^2*Id
    P[N*(N-1)+1:N^2, N*(N-1)+1:N^2] += -hi^2*Id
    #P[1, :] .= 0
    #P[1, 1] = 1*hi^2
    #P[N^2, :] .= 0
    #P[N^2, 1] = 1*hi^2
    return P
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

    ω[N+2, :] = 2/h^2 .* (rn - rn1) .- 2*u_top/h
    ω[1, :] = 2/h^2 .* (r3 - r2) .+ 2*u_bottom/h

    ω[:, 1] = 2/h^2 .* (c3 - c2) .- 2*v_left/h
    ω[:, N+2] = 2/h^2 .* (cn - cn1) .+ 2*v_right/h
   
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

    for i in 3:N
        for j in 3:N
            # higher indices currently correspond to top
            # u = dψ/dy
            u_ji = 1/(2*h)*(ψ[j+1, i] - ψ[j-1, i])
            # v = -dψ/dx
            v_ji = -1/(2*h)*(ψ[j, i+1] - ψ[j, i-1])
            ωx = (u_ji >= 0) ? 
                    (1/(3*h)*(-ω[j, i+2] + 3*ω[j, i+1] - 3*ω[j, i] + ω[j, i-1])) :
                    (1/(3*h)*(-ω[j, i+1] + 3*ω[j, i] - 3*ω[j, i-1] + ω[j, i-2]))
            ωy = (v_ji >= 0) ?
                    (1/(3*h)*(-ω[j+2, i] + 3*ω[j+1, i] - 3*ω[j, i] + ω[j-1, i])) :
                    (1/(3*h)*(-ω[j+1, i] + 3*ω[j, i] - 3*ω[j-1, i] + ω[j-2, i]))

            ω_next[j, i] = ω[j, i] + 
                        Δt*( -1/(2*h)*u_ji*(ω[j, i+1] - ω[j, i-1])
                            -1/(2*h)*v_ji*(ω[j+1, i] - ω[j-1, i])
                            + 1/(Re*h^2)*(ω[j, i+1] - ω[j, i-1] - 4*ω[j, i] + ω[j+1, i] + ω[j-1, i])
                            - 1/2*u_ji*ωx
                            - 1/2*v_ji*ωy
                            )
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
    ω = zeros(N+2, N+2)
    ωtmp = zeros(N+2, N+2)
    ψ = zeros(N+2, N+2)
    P = genPoissonMtx(opts)
    ϵ = 1e-10

    max_iter = 10

    for i in max_iter
        updateVorticityWallBC!(ω, ψ, opts, opts_BC)
        ω  = updateVorticity(ω, ωtmp, ψ, opts)
        ω_inner = @view ω[2:end-1, 2:end-1]
        ω_vec = reshape(ω_inner, N*N, 1)
        poisson_RHS = -ω_vec
        ψ_inner = @view ψ[2:end-1, 2:end-1]
        ψ_vec = reshape(ψ_inner, N*N, 1)
        ψ, num_iter = CG_std(P, poisson_RHS, ψ_vec, ϵ)
    end

end

n = 64
h = 1/N
dt = 0.001
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
            "Re"=>150
            )
runVorticityStream(opts, BC_opts)