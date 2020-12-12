
# NOTE! the signs may be incorrect!!
# double check
function predictorSIMPLE(u, v, p, u_new, v_new, opts)
    Δt = opts["dt"]
    Δx, Δy = opts["dx"], opts["dy"]
    Re_inv = opts["Rei"]

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
end

function velCorrectorSIMPLE(u_corrector, v_corrector, p_corrector, opts)
    Δt = opts["dt"]
    Δx, Δy = opts["dx"], opts["dy"]
    Δu = 1 / (Δt*Δx) * (Δp_right - Δp)
    Δv = 1 / (Δt*Δy) * (Δp_down - Δp)
    return Δu, Δv
end