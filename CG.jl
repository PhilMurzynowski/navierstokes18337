using LinearAlgebra, Printf

"""
Included methods in CG.jl
Conjugate Gradient                      CG(A, b, x_guess, ϵ, max_iter=1e3)
Preconditioned CG                       PCG(A, Minv, b, x_guess, ϵ, max_iter=1e3)
CG specifcally for Poisson              CG_Poisson(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3)
    (no A matrix passed in)
Incomplete Cholesky CG                  ICCG(A, U, b, x_guess, ϵ, max_iter=1e3)
    (fastest)
    (may be updated to use preallocated tmp arrays as well)
"""

"""
CG, Conjugate Gradients implementation
Symmetric Matrix A must be constructed and passed in.
Warning:
    Do not use tolerance very close to machine error
    as this funciton does not specially account for roundoff error.
"""
function CG(A, b, x_guess, ϵ, max_iter=1e3)

    residual = b - A*x_guess
    search_direction = residual     # use residuals as conjugate search directions
    rTr = dot(residual, residual)   # save in variable to avoid repeat calculations
    rTr_next = nothing              # will need 2 vars for rTr
    x = x_guess
    iter = 0

    while rTr > ϵ^2 && iter < max_iter
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

"""
Same as CG_std except also making use of preconditioning
Minv is the preconditioner, passsed in as a matrix.
Now requires an extra matrix vector product for preconditioning
    each iteration.
Note: 
    Consider modifying the function passing in preallocated memory.
"""
function PCG(A, Minv, b, x_guess, ϵ, max_iter=1e3)

    residual = b - A*x_guess
    presidual = Minv*residual
    search_direction = presidual     # use preconditioned residuals as conjugate search directions
    rTpr = dot(residual, presidual)  # save in variable to avoid repeat calculations
    rTpr_next = nothing              # will need 2 vars for rTr
    x = x_guess
    iter = 0

    while rTpr > ϵ^2 && iter < max_iter
        iter += 1
        mvp = A*search_direction
        step_size = rTpr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        presidual = Minv*residual
        rTpr_next = dot(residual, presidual)
        gsc = rTpr_next / rTpr        # gram-schmidt elimination
        rTpr = rTpr_next
        # update x
        x += step_size.*search_direction
        search_direction = presidual + gsc.*search_direction
    end

    return x, iter # return solution and number of iterations for analysis

end

"""
ICCG : Incomplte Cholesky Conjugate Gradient
Very similar to PCG except using Incomplete Cholesky
for preconditioning, so can make more optimizations,
i.e.
    Using faster solves for banded, sparse triangular matrices instead of calculating
    Minv which may be dense.

U is the Upper triangular of the incomplete Cholesky factorization

UPDATE TO USE PREALLOCATED ARRAYS
"""
function ICCG(A, U, b, x_guess, ϵ, max_iter=1e3)

    residual = b - A*x_guess
    # two triangular solves
    # Faster than multiplying by Minv
    presidual = U' \ residual
    presidual = U \ presidual           # preconditioned residual

    search_direction = presidual        # use preconditioned residuals as conjugate search directions
    rTpr = dot(residual, presidual)     # save in variable to avoid repeat calculations
    rTpr_next = nothing                 # will need 2 vars for rTr
    x = x_guess
    iter = 0

    # Warning: since using rTpr, not rTr, approximate
    while rTpr > ϵ^2 && iter < max_iter
        #println(iter)
        iter += 1
        mvp = A*search_direction
        step_size = rTpr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        presidual = U' \ residual
        presidual = U \ presidual
        rTpr_next = dot(residual, presidual)
        gsc = rTpr_next / rTpr        # gram-schmidt elimination
        rTpr = rTpr_next
        # update x
        x += step_size.*search_direction
        search_direction = presidual + gsc.*search_direction
    end

    return x, iter # return solution and number of iterations for analysis

end

# test CG_std
function test_GC(A=nothing, b=nothing, ϵ=1e-9)

    if A === nothing || b === nothing
        # very likely to be full rank if random
        # not specifically testing sparse solving ability at the moment
        A = rand(100, 100)
        A = A' + A # make symmetric
        b = rand(100)
    end

    x_direct = A \ b
    x_CG, num_iter = CG(A, b, vec(zeros(length(b))), ϵ)
    @printf "number iterations: %d for epsilon: %.16f\n" num_iter ϵ
    #display(x_direct)
    #display(x_direct)
    @assert isapprox(x_direct, x_CG, atol=ϵ)
end

#test_GC_std()

"""
CG gradient method specialized for 2D Poisson equation
with Dirchlet Boundary conditions on all 4 sides.

Warning: Currently buggy, has also been replaced by ICCG.

Main advantage: Do not need to generate and pass in Poisson matrix.
Main disadvantages: cannot easily take advantage of fast matrix mulitplies in Julia
    for sparse, banded matrices, and cannot easily apply a general class of preconditioners,
    not scalable to write out all multiplications in loops.

Due to the disadvantages above, CG, PCG, and ICCG preferred if using special sparse or banded matrices.

Notes: Everything is kept and expected to be passed in a 2D Array, makes for more
        readable indexing.
    Uses preallocated memory.
Optimization Notes:
    A new tmp variable can be used so that x does not have to be sliced every iteration.
    Only slice and update x just before returning it
        x[2:end-1, 2:end-1] .+= step_size*search_direction
"""
function CG_Poisson(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3)
    N = opts["N"]
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
    # only one slice, so don't slice during iteration
    #x = x_guess[2:end-1, 2:end-1]
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

    while rTr > ϵ^2
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

        # will += of a slice create a new array... check
        x[2:end-1, 2:end-1] .+= step_size*search_direction
        search_direction = residual + gsc.*search_direction
    end

    return x, iter

    #=
    TODO: optimization
        # using broadcast operation instead of loop to encourage vectorizing
        x .+= step_size*search_direction
        search_direction = residual + gsc.*search_direction
    end

    # reuse x_guess
    # Note: chance that this is buggy! double check
    x_guess[2:end-1, 2:end-1] = x

    return x_guess, iter
    =#

end

"""
Warning: This is actually totally useless precisely because all the diagonal entries
are the same, reasoned after, now kept as proof of concept.

Warning: Currently buggy, has also been replaced by ICCG.

Same as CG_Poisson except also optimized to use a simple diagonal proconditioner.
Furthermore, since the 2D Poisson Matrix has the same constant value on the diagonals,
it is especially easy to implement using broadcast operations.

Optimization Notes:
    Apply M^-1 to the initual residual once instead of both A and b,
    and keep applying to residual.
"""
function PCG_Poisson_diagonal(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3, tmp4)
    N = opts["N"]
    h = opts["h"]
    ih_sq = 1/h^2

    # know diagonal values from Poisson mtx structure
    idiag = -h^2/4
    idiag = 1.0

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

    # approximate, not exactly rTr
    while rTpr > ϵ^2
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
        #@printf "iter: %d, rTr: %3.10f\n" iter rTpr
    end

    return x, iter

end