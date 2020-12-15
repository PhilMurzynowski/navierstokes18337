using LinearAlgebra, Printf

# CG, Conjugate Gradients implementation

# standard method involving multiplying matrices
# matrix A must be constructed and passed in
# warning do not use tolerance very close to machine error
#   as this funciton does not specially account for roundoff error
function CG_std(A, b, x_guess, ϵ)

    residual = b - A*x_guess
    search_direction = residual     # use residuals as conjugate search directions
    rTr = dot(residual, residual)   # save in variable to avoid repeat calculations
    rTr_next = nothing              # will need 2 vars for rTr
    x = x_guess
    iter = 0

    # using l1 norm so don't have to square tiny ϵ
    while sum(abs.(residual)) > ϵ
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
    @printf "number iterations: %d for epsilon: %f\n" num_iter ϵ
    #display(x_direct)
    #display(x_direct)
    @assert isapprox(x_direct, x_CG, atol=ϵ)
end

test_GC_std()