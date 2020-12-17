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
function CG_Poisson(b, p_guess, ϵ, opts)
    N = opts["N"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    b = reshape(b, N, N)
    # compute initial residual
    # can pass in in the future
    residual = zeros(N, N)
    for i in 2:N-1
        #i_inner = i-1  # for indexing convenience
        for j in 2:N-1
            #j_inner = j - 1
            #residual[j, i] = b[j_inner, i_inner] - 1/Δx^2*(p_guess[j, i+1] - 2*p_guess[j, i] + p_guess[j, i-1])
            residual[j, i] = b[j, i] - 1/Δx^2*(p_guess[j, i+1] - 2*p_guess[j, i] + p_guess[j, i-1])
            residual[j, i] -= 1/Δy^2*(p_guess[j+1, i] - 2*p_guess[j, i] + p_guess[j-1, i])
        end
    end
    # any bugs with copying here?
    search_direction = residual
    rTr = dot(residual, residual)
    rTr_next = nothing
    pressure = p_guess
    iter = 0
    # pass in in the future to reuse
    mvp = zeros(N, N) # init storage for matrix vector product
    
    # using l1 norm so don't have to square tiny ϵ
    while norm(residual, 1) > ϵ
        iter += 1
        # single matrix vector product optimization
        for i in 2:N-1
            for j in 2:N-1
                mvp[j, i] = 1/Δx^2*(search_direction[j, i+1] - 2*search_direction[j, i] + search_direction[j, i-1])
                mvp[j, i] += 1/Δy^2*(search_direction[j+1, i] - 2*search_direction[j, i] + search_direction[j-1, i])
            end
        end
        # calculate step size
        step_size = rTr / dot(search_direction, mvp)
        residual -= step_size.*mvp
        rTr_next = dot(residual, residual)
        gsc = rTr_next / rTr
        rTr = rTr_next
        # update pressure solution
        # One way to update pressure
        # TODO: However with with a better iteration pattern hopefully avoid cachemisses
        pressure += step_size.*search_direction
        
        search_direction = residual + gsc.*search_direction
    end

    return pressure, iter

end