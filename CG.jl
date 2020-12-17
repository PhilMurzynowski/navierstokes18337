using LinearAlgebra, Printf

# CG, Conjugate Gradients implementation

# standard method involving multiplying matrices
# matrix A must be constructed and passed in
# warning do not use tolerance very close to machine error
#   as this funciton does not specially account for roundoff error
function CG_std(A, b, x_guess, ϵ, max_iter=1e3)

    residual = b - A*x_guess
    search_direction = residual     # use residuals as conjugate search directions
    rTr = dot(residual, residual)   # save in variable to avoid repeat calculations
    rTr_next = nothing              # will need 2 vars for rTr
    x = x_guess
    iter = 0

    # using l1 norm so don't have to square tiny ϵ
    #while sum(abs.(residual)) > ϵ
    while norm(residual, 1) > ϵ && iter < max_iter
        #println(iter)
        iter += 1
        # single matrix vector product optimization
        mvp = A*search_direction
        step_size = rTr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        rTr_next = dot(residual, residual)
        gsc = rTr_next / rTr        # gram-schmidt elimination
        rTr = rTr_next
        # update x
        x += step_size.*search_direction
        search_direction = residual + gsc.*search_direction
    end

    return x, iter # return solution and number of iterations for analysis

end

# test CG_std
function test_GC_std(A=nothing, b=nothing, ϵ=1e-9)

    if A === nothing || b === nothing
        # very likely to be full rank if random
        # not specifically testing sparse solving ability at the moment
        A = rand(100, 100)
        A = A' + A # make symmetric
        b = rand(100)
    end

    x_direct = A \ b
    x_CG, num_iter = CG_std(A, b, vec(zeros(length(b))), ϵ)
    @printf "number iterations: %d for epsilon: %.16f\n" num_iter ϵ
    #display(x_direct)
    #display(x_direct)
    @assert isapprox(x_direct, x_CG, atol=ϵ)
end

#test_GC_std()

"""
CG gradient method specialized for 2D Poisson equation
with Dirchlet Boundary conditions on all 4 sides.

Do not need to generate Poisson matrix.

Note: everything is kept and expected to be passed in a 2D Array, makes for more
readable indexing
"""
function CG_Poisson(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3)
    N = opts["N"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    h = opts["h"]
    ih_sq = 1/h^2

    # compute initial residual
    # residual should be N-2, N-2
    # can keep it slightly smaller as only iterating over inner points
    # leads to some slight loop unrolling and complexity later
    # does not matter if garbage in tmp1,
    # residual will be fully overwritten before use 
    residual = tmp1

    @inbounds for i in 2:N-1
        i_inner = i-1  # for indexing convenience
        @inbounds for j in 2:N-1
            j_inner = j - 1
            residual[j_inner, i_inner] = b[j, i] - ih_sq*(x_guess[j, i+1] - 2*x_guess[j, i] + x_guess[j, i-1])
            residual[j_inner, i_inner] -= ih_sq*(x_guess[j+1, i] - 2*x_guess[j, i] + x_guess[j-1, i])
        end
    end
    # same idea with search_direction, reuse array
    # have to copy residual once initially
    search_direction = tmp2
    search_direction = copy(residual)
    rTr = dot(residual, residual)
    rTr_next = nothing
    x = x_guess
    iter = 0
    # storage for matrix vector product
    # passed memory for reuse
    # does not matter if garbage in tmp1,
    # mvp will be fully overwritten before use
    mvp = tmp3

    
    @inline bottomleft(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j+1, i])
    @inline left(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j+1, i] + sd[j-1, i])
    @inline topleft(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j-1, i])
    @inline inner(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j+1, i] + sd[j-1, i])
    @inline bottomright(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j+1, i])
    @inline right(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j+1, i] + sd[j-1, i])
    @inline topright(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j-1, i])
    @inline top(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j-1, i])
    @inline bottom(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j+1, i])

    # using l1 norm so don't have to square tiny ϵ
    while norm(residual, 1) > ϵ
        iter += 1

        # left column
        mvp[1, 1] = bottomleft(search_direction, 1, 1)
        @inbounds for j in 2:N-3 
            mvp[j, 1] = left(search_direction, j, 1)
        end
        mvp[N-2, 1] = topleft(search_direction, N-2, 1)
        # end of left column

        # inner points
        @inbounds for i in 2:N-3
            j = 1
            mvp[j, i] = bottom(search_direction, j, i)
            @inbounds for j in 2:N-3
                mvp[j, i] = inner(search_direction, j, i)
            end
            j = N-2
            mvp[j, i] = top(search_direction, j, i)
        end
        # end of inner points

        # right column
        mvp[1, N-2] = bottomright(search_direction, 1, N-2)
        @inbounds for j in 2:N-3
            mvp[j, N-2] = right(search_direction, j, N-2)
        end
        mvp[N-2, N-2] = topright(search_direction, N-2, N-2)
        # end of right column

        # calculate step size
        step_size = rTr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        rTr_next = dot(residual, residual)
        gsc = rTr_next / rTr
        rTr = rTr_next
        # update x solution
        # using broadcast operation instead of loop to encourage vectorizing
        # will += of a slice create a new array... check
        x[2:end-1, 2:end-1] .+= step_size*search_direction
        search_direction = residual + gsc.*search_direction
    end

    return x, iter

end

"""
Same as CG_Poisson except also optimized to use a simple diagonal proconditioner.
Furthermore, since the 2D Poisson Matrix has the same constant value on the diagonals,
it is especially easy to implement using broadcast operations.

Apply M^-1 to the initual residual once instead of both A and b,
and keep applying to residual.
"""
function PCG_Poisson_diag(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3, tmp4)
    N = opts["N"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    h = opts["h"]
    ih_sq = 1/h^2

    # know diagonal values from Poisson mtx structure
    idiag = -h^2/4

    # compute initial residual
    # residual should be N-2, N-2
    # can keep it slightly smaller as only iterating over inner points
    # leads to some slight loop unrolling and complexity later
    # does not matter if garbage in tmp1,
    # residual will be fully overwritten before use 
    residual = tmp1

    @inbounds for i in 2:N-1
        i_inner = i-1  # for indexing convenience
        @inbounds for j in 2:N-1
            j_inner = j - 1
            residual[j_inner, i_inner] = b[j, i] - ih_sq*(x_guess[j, i+1] - 2*x_guess[j, i] + x_guess[j, i-1])
            residual[j_inner, i_inner] -= ih_sq*(x_guess[j+1, i] - 2*x_guess[j, i] + x_guess[j-1, i])
        end
    end
    # same idea, reuse array
    presidual = tmp4    # using preconditioner
    presidual = residual .* idiag


    # have to copy residual once initially
    search_direction = tmp2
    search_direction = copy(presidual)
    # instead of doing the first commented out operation below
    # can dot product with itself so that it hopefully
    # uses half as much memory, and mulitply once after
    #rTpr = dot(residual, presidual)
    rTpr = dot(residual, residual) * idiag
    rTpr_next = nothing
    x = x_guess
    iter = 0
    # storage for matrix vector product
    # passed memory for reuse
    # does not matter if garbage in tmp1,
    # mvp will be fully overwritten before use
    mvp = tmp3

    
    @inline bottomleft(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j+1, i])
    @inline left(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j+1, i] + sd[j-1, i])
    @inline topleft(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j-1, i])
    @inline inner(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j+1, i] + sd[j-1, i])
    @inline bottomright(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j+1, i])
    @inline right(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j+1, i] + sd[j-1, i])
    @inline topright(sd, j, i) = ih_sq*(-4*sd[j, i] + sd[j, i-1] + sd[j-1, i])
    @inline top(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j-1, i])
    @inline bottom(sd, j, i) = ih_sq*(sd[j, i+1] - 4*sd[j, i] + sd[j, i-1] + sd[j+1, i])

    # using l1 norm so don't have to square tiny ϵ
    while norm(residual, 1) > ϵ
        iter += 1

        # left column
        mvp[1, 1] = bottomleft(search_direction, 1, 1)
        @inbounds for j in 2:N-3 
            mvp[j, 1] = left(search_direction, j, 1)
        end
        mvp[N-2, 1] = topleft(search_direction, N-2, 1)
        # end of left column

        # inner points
        @inbounds for i in 2:N-3
            j = 1
            mvp[j, i] = bottom(search_direction, j, i)
            @inbounds for j in 2:N-3
                mvp[j, i] = inner(search_direction, j, i)
            end
            j = N-2
            mvp[j, i] = top(search_direction, j, i)
        end
        # end of inner points

        # right column
        mvp[1, N-2] = bottomright(search_direction, 1, N-2)
        @inbounds for j in 2:N-3
            mvp[j, N-2] = right(search_direction, j, N-2)
        end
        mvp[N-2, N-2] = topright(search_direction, N-2, N-2)
        # end of right column

        # calculate step size
        step_size = rTpr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        presidual = residual .* idiag
        rTpr_next = dot(residual, residual) * idiag
        gsc = rTpr_next / rTpr
        rTpr = rTpr_next
        # update x solution
        # using broadcast operation instead of loop to encourage vectorizing
        # will += of a slice create a new array... check
        x[2:end-1, 2:end-1] .+= step_size*search_direction
        search_direction = presidual + gsc.*search_direction
    end

    return x, iter

end