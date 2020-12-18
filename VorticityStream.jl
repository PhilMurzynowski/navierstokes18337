using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap, contour
using Printf

# example included at end of file

include("CG.jl")
include("Matrices.jl")


# could this be fused to updateVoriticity for better cache use, repeating this question everywhere ha
function vs_updateVorticityWallBC!(ω, ψ, opts, opts_BC)

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
        ω[j, N] = -2/h^2*ψ[j, N-1] + 2/h*v_right
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
"""
Update vorticity one timestep
Note: use normal edition for more performance
      left for benchmarking convenience against standard version
"""
function vs_updateVorticity_branching(ω, ω_next, ψ, opts)
    Δt = opts["dt"]
    h = opts["h"]   # just use h since know using dx and dy equal
    Re = opts["Re"]
    # large j correspond to higher point in grid
    for i in 2:N-1
        i_boundary = (i == 2 || i == N-1) ? true : false
        for j in 2:N-1
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
            if i_boundary || j == 2 || j == N-1
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
        end
    end
    return ω_next
end


"""
Update vorticity one timestep
Notes:
    @inbounds may likely be overkill given recent Julia updates, included for assurance
"""
function vs_updateVorticity(ω, ω_next, ψ, opts)
    Δt = opts["dt"]
    h = opts["h"]   # just use h since know using dx and dy equal
    Re = opts["Re"]

    # pulled some constants with division out
    # might be taken care of by compiler anyway
    ih2 = 1/(2*h)
    ih3 = 1/(3*h)
    ih_sq = 1/(h^2)
    iRe = 1/Re

    # inline functions as they're just a few expressions
    # standard update, used on boundaries
    @inline function update_std(j, i)
        u_ji = ih2*(ψ[j+1, i] - ψ[j-1, i]) # u = dψ/dy
        v_ji = -ih2*(ψ[j, i+1] - ψ[j, i-1]) # v = -dψ/dx
        ω_next[j, i] = ω[j, i] + 
                    Δt*( -ih2*u_ji*(ω[j, i+1] - ω[j, i-1])
                        -ih2*v_ji*(ω[j+1, i] - ω[j-1, i])
                        + ih_sq/Re*(ω[j, i+1] - ω[j, i-1] - 4*ω[j, i] + ω[j+1, i] + ω[j-1, i]))
    end
    # if not on boundary can use 2nd order upwind for more accuracy
    @inline function update_order2_upwind(j, i)
        u_ji = ih2*(ψ[j+1, i] - ψ[j-1, i]) # u = dψ/dy
        v_ji = -ih2*(ψ[j, i+1] - ψ[j, i-1]) # v = -dψ/dx

        # use ternary assignment instead of full branching, should be faster
        ωx = (u_ji >= 0) ? 
                ih3*(-ω[j, i+2] + 3*ω[j, i+1] - 3*ω[j, i] + ω[j, i-1]) :
                ih3*(-ω[j, i+1] + 3*ω[j, i] - 3*ω[j, i-1] + ω[j, i-2])
        ωy = (v_ji >= 0) ?
                ih3*(-ω[j+2, i] + 3*ω[j+1, i] - 3*ω[j, i] + ω[j-1, i]) :
                ih3*(-ω[j+1, i] + 3*ω[j, i] - 3*ω[j-1, i] + ω[j-2, i])

        ω_next[j, i] = ω[j, i] + 
                    Δt*( -ih2*u_ji*(ω[j, i+1] - ω[j, i-1])
                        -ih2*v_ji*(ω[j+1, i] - ω[j-1, i])
                        + ih_sq*iRe*(ω[j, i+1] - ω[j, i-1] - 4*ω[j, i] + ω[j+1, i] + ω[j-1, i]))
        ω_next[j, i] += Δt*(
                        - 1/2*u_ji*ωx
                        - 1/2*v_ji*ωy
                        )
    end

    # left boundary
    i = 2
    @inbounds for j in 2:N-1
        update_std(j, i)
    end
    # end of left boundary code
    @inbounds for i in 3:N-2
        # bottom boundary
        j = 2
        update_std(j, i)
        # end of bottom boundary
        @inbounds for j in 3:N-2
            update_order2_upwind(j, i)
        end
        # bottom boundary
        j = N-1
        update_std(j, i)
        # end of bottom boundary
    end
    # right boundary
    i = N-1
    @inbounds for j in 2:N-1
        update_std(j, i)
    end
    # end of right boundary code

    return ω_next
end

# solve poisson with CG

# is computing velocity necessary?
# believe can eliminate storing velocity


function run_vorticitystream_simulation(opts, opts_BC)
    N = opts["N"]
    ϵ = opts["ϵ"]
    timesteps = opts["timesteps"]

    # init vorticity, streamfunction
    ω = zeros(N, N)
    ωtmp = zeros(N, N)
    ψ = zeros(N, N)

    # arrays for memory reuse, preallocating
    tmp1 = zeros(N-2, N-2)
    tmp2 = zeros(N-2, N-2)
    tmp3 = zeros(N-2, N-2)
    tmp4 = zeros(N-2, N-2)
    # use if not using CG_Poisson
    use_special = false
    if use_special == false
        P = genPoissonMtx(N-2, opts)
        # preconditioner
        # sanity check with identity
        Minv = Diagonal(ones(size(P)))
        # Incomplete Cholesky
        chol = cholesky(P)
        U = chol.U
        U[P .== 0.0] .= 0
        #Minv = inv(U'*U)
    end

    for i in 1:timesteps
        vs_updateVorticityWallBC!(ω, ψ, opts, opts_BC)
        ω  = vs_updateVorticity(ω, ωtmp, ψ, opts)
        # perhaps keep vorticity negated if going to be flipping sign
        # or do the negative in the CG_Poisson solver

        # optimize away if statment later when not comparing
        if use_special
            # could modify this further so don't have to negate vorticity
            # when passing it in, simply keep negative vorticity as a variable
            # but negating not expensive, and could lead to bugs down the line
            #ψ_copy = copy(ψ)
            ψ, num_iter1 = CG_Poisson(-ω, ψ, ϵ, opts, tmp1, tmp2, tmp3)
            #ψ, num_iter2 = PCG_Poisson_diag(-ω, ψ, ϵ, opts, tmp1, tmp2, tmp3, tmp4)
            #@printf "CG: %d, PCG_diag %d\n" num_iter1 num_iter2
            #@printf "PCG_diag %d\n" num_iter2
            @printf "CG %d\n" num_iter1
        else
            ω_inner = @view ω[2:end-1, 2:end-1]
            ω_vec = reshape(ω_inner, (N-2)*(N-2), 1)
            # double check if negative or not
            poisson_RHS = ω_vec
            #poisson_RHS = -ω_vec
            ψ_inner = @view ψ[2:end-1, 2:end-1]
            ψ_vec = vec(ψ_inner)

            #ψ_vec_copy = copy(ψ_vec)
            #poisson_RHS_copy = copy(poisson_RHS)
            #ψ_vec_copy, num_iter1 = CG_std(P, poisson_RHS_copy, ψ_vec_copy, ϵ, N^2)
            #ψ_vec, num_iter2 = PCG_std(P, Minv, poisson_RHS, ψ_vec, ϵ, N^2)
            ψ_vec, num_iter2 = ICCG(P, U, poisson_RHS, ψ_vec, ϵ, N^2)
            #@printf "CG: %d, PCG %d\n" num_iter1 num_iter2
            #ψ_vec, num_iter = ICCG_std(P, L, poisson_RHS, ψ_vec, ϵ, N^2)
            #ψ_vec = P \ poisson_RHS
            ψ[2:end-1, 2:end-1] = reshape(ψ_vec, N-2, N-2)
        end
    end

    return ω, ψ
end

function plot_ωψ(ω, ψ)
    h = opts["h"]
    N = opts["N"]
    xs = 0.0:h:h*(N-1)
    ys = 0.0:h:h*(N-1)
    c = contour(xs, ys, ψ')
end

#=
# example

N = 32
h = 1/N
dt = 0.001
BC_opts = Dict("u_top"=>1.0,
               "u_bottom"=>0.0,
               "v_left"=>0.0,
               "v_right"=>1.0)
opts = Dict("N"=>N,
            "Nx"=>N,
            "Ny"=>N,
            "dx"=>h,
            "dy"=>h,
            "h"=>h,
            "dt"=>dt,
            "Re"=>10.0,
            "ϵ"=>1e-4,         # tolerance parameter for CG
            "timesteps"=>100
            )
@time ω, ψ = run_vorticitystream_simulation(opts, BC_opts)
plot_ωψ(ω, ψ)

=#