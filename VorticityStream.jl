


# init

# compute BC

# simulate vorticity forward a timestep
# starting with simple FTCS (forward time centered space) first order method
#   2nd order upwind for convection
#   check if worth over 1st order, not much harder to implement
# can switch to Dormand-Prince
# should I be storing u v in arrays
# does second order differencing have more potential for cache misses? consequence
function updateVorticity(ω, u, v, ω_next, opts)
    Δt = opts["dt"]
    h = opts["h"]   # just use h since know using dx and dy equal
    Re = opts["Re"]

    for i in 1:N
        for j in 1:N
            ωx = (u[j, i] >= 0) ? 
                    (1/(3*h)*(-ω[j, i+2] + 3*w[j, i+1] - 3*ω[j, i] + ω[j, i-1])) :
                    (1/(3*h)*(-ω[j, i+1] + 3*w[j, i] - 3*ω[j, i-1] + ω[j, i-2]))
            ωy = (v[j, i] >= 0) ?
                    (1/(3*h)*(-ω[j+2, i] + 3*w[j+1, i] - 3*ω[j, i] + ω[j-1, i])) :
                    (1/(3*h)*(-ω[j+1, i] + 3*w[j, i] - 3*ω[j-1, i] + ω[j-2, i]))

            ω_next[j, i] = ω[j, i] + 
                        Δt*( -1/(2*h)*u[i, j]*(ω[j, i+1] - ω[j, i-1])
                            -1/(2*h)*v[i, j]*(ω[j+1, i] - ω[j-1, i])
                            + 1/(Re*h^2)*(ω[j, i+1] - ω[j, i-1] - 4*ω[j, i] + ω[j+1, i] + ω[j-1, i])
                            - 1/2*u[j, i]*ωx
                            - 1/2*v[j, i]*ωy
                            )
        end
    end
    return
end

# could this be fused to updateVoriticity for better cache use, repeating this question everywhere ha
function updateVorticityWallBC(ω, u, v, ψ, opts, BC_opts)
    h = opts["h"]

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

    ω[N+2, :] = 2/h^2 .* (rn - rn1) - 2*u_top/h
    ω[1, :] = 2/h^2 .* (r3 - r2) - 2*u_bottom/h
    w[i,1] = (-4.0*s[i,2]+0.5*s[i,3])/(dy*dy)

end

# solve poisson with CG

# is computing velocity necessary?