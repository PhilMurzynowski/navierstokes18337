


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

    for i in 2:N-1
    for j in 2:N-1
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
    return
end

# solve poisson with CG

# is computing velocity necessary?