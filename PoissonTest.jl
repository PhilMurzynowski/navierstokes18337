using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap
using BenchmarkTools

include("Matrices.jl")
include("CG.jl")
include("gridSolver.jl")

"""
File dedicated to testing and benchmarking different approaches to solving Poisson Eq.
Testing approach:
    1. A solution, u, is generated (e.g. an example pressure field)
    2. Source terms are determined by multiplying u by laplacian
    3. Validate solver by comparing to starting u after solving Poisson Eq given source terms
Note: A nice update would be to use surface plots instead later! would be cooler!
"""

"""
Given a u, return the appropriate source vector
"""
function determine_source(u, P)
    source = P*vec(u)
    return source
end

"""
Given u, create heatmap
"""
function get_u_scene(xs, ys, u)
    # likely don't need this reverse tranpose shenangians but ok for now
    hm = heatmap(xs, ys, reverse(u', dims=2), colormap=:berlin, show_axis=false)
    cl = colorlegend(hm[end], raw = true, camera = campixel!)
    return hm, cl
end

"""
Given u show that all solvers output the correct solution.
"""
function visual_test(u, xs, ys, opts)
    N = opts["N"]
    h = opts["h"]
    ϵ = opts["ϵ"]
    hm_actual, cl_actual = get_u_scene(xs, ys, u)

    # UPDATE genPoissonMtx
    P = genPoissonMtx(N, h)
    # Incomplete Cholesky
    # UDPATE make all of these banded!
    chol = cholesky(P)
    U = chol.U
    U[P .== 0.0] .= 0
    Minv = inv(U'*U)

    source = determine_source(u, P)

    # zero out every time to be sure
    x_guess = zeros(length(source))
    CG_sol, iter = CG(P, source, x_guess, ϵ)
    hm_CG, cl_CG = get_u_scene(xs, ys, reshape(CG_sol, N, N))

    x_guess = zeros(length(source))
    PCG_sol, iter = PCG(P, Minv, source, x_guess, ϵ)
    hm_PCG, cl_PCG = get_u_scene(xs, ys, reshape(PCG_sol, N, N))
    
    x_guess = zeros(length(source))
    ICCG_sol, iter = ICCG(P, U, source, x_guess, ϵ)
    hm_ICCG, cl_ICCG = get_u_scene(xs, ys, reshape(ICCG_sol, N, N))

    #=
    x_guess = zeros(N, N)
    tmp1 = zeros(N-2, N-2)
    tmp2 = copy(tmp1)
    tmp3 = copy(tmp1)
    tmp4 = copy(tmp1)
    # note, meant to use CG_Poisson here, but broken for some reason and PCG_Poisson_diag is equivalent
    CGPoisson_sol, iter = PCG_Poisson_diagonal(reshape(source, N, N), x_guess, ϵ, opts, tmp1, tmp2, tmp3, tmp4)
    hm_CGPoisson, cl_CGPoisson = get_u_scene(xs, ys, CGPoisson_sol)
    =#

    x_guess = zeros(N, N)
    tmp = copy(x_guess)
    Grid_sol, iter = grid_poisson_solve(x_guess, tmp, reshape(source, N, N)[2:end-1, 2:end-1], opts, 1e3)
    hm_grid, cl_grid = get_u_scene(xs, ys, Grid_sol)

    # hacky way to arrange things
    parent = Scene(resolution= (800, 400))
    full_scene = hbox(
                    vbox(hm_grid, cl_grid), #note shows up on bottom when plotting
                    vbox(hm_actual, cl_actual),
                    vbox(hm_CG, cl_CG), 
                    vbox(hm_PCG, cl_PCG), 
                    #vbox(hm_CGPoisson, cl_CGPoisson), 
                    vbox(hm_ICCG, cl_ICCG), 
                    parent=parent)

    # separate plot if buggy
    #=
    parent = Scene(resolution= (300, 300))
    grid_scene = hbox(
                    vbox(hm_grid, cl_grid), 
                    parent=parent)
    =#
    #println(iter)
    display(full_scene)
    return
end

"""
Given a source vector compare timings of different solution methods.
"""
function timing_test(u, opts)

    N = opts["N"]
    h = opts["h"]
    #ϵ = opts["ϵ"]

    Ptmp = genPoissonMtx(N, h)
    # Incomplete Cholesky
    chol = cholesky(Ptmp)
    Utmp = chol.U
    Utmp[Ptmp .== 0.0] .= 0
    Minv = inv(Utmp'*Utmp)
    # convert both to sparse
    # doing this workaround because of bugs with setindex! and setting to 0
    # labeling U as UpperTriangular seems to help a smidgen
    P = sparse(Ptmp)
    U = UpperTriangular(sparse(Utmp))

    source = determine_source(u, P)

    datapts = 20
    ϵs = 10 .^ range(-4, -10, length=datapts)
    CG_time = zeros(datapts)
    PCG_time = zeros(datapts)
    ICCG_time = zeros(datapts)
    grid_time = zeros(datapts)

    CG_iter = zeros(datapts)
    PCG_iter = zeros(datapts)
    ICCG_iter = zeros(datapts)
    grid_iter = zeros(datapts)

    i = 0
    for j in 1:datapts+1

        i+= 1
        ϵ = ϵs[i]

        # randomizing x_guess to prevent optimization
        x_guess = rand(-1e-3:1e3, length(source))

        CG_time[i] = @elapsed CG_sol, CG_iter[i] = CG(P, source, x_guess, ϵ)

        PCG_time[i] = @elapsed PCG_sol, PCG_iter[i] = PCG(P, Minv, source, x_guess, ϵ)
        
        ICCG_time[i] = @elapsed ICCG_sol, ICCG_iter[i] = ICCG(P, U, source, x_guess, ϵ)

        # note, meant to use CG_Poisson here, but broken for some reason

        tmp = reshape(copy(x_guess), N, N)
        grid_time[i] = @elapsed grid_sol, grid_iter[i] = grid_poisson_solve(reshape(x_guess, N, N), tmp, reshape(source, N, N)[2:end-1, 2:end-1], opts, 2e3)

        @printf "Size (N): %d, ϵ: %.10f\n" N ϵ
        @printf "Times:\n Grid: %f, CG: %f, PCG: %f, ICCG: %f\n" grid_time[i] CG_time[i] PCG_time[i] ICCG_time[i]
        @printf "Iterations:\n Grid: %d, CG: %d, PCG: %d, ICCG: %d\n" grid_iter[i] CG_iter[i] PCG_iter[i] ICCG_iter[i]

        # from 2nd run on will be faster to to JIT compile
        # so overwrite first
        if j == 1
            i -= 1
        end
    end

    timings = Scene()
    eps = reverse(ϵs)
    plot!(eps, grid_time)
    plot!(eps, CG_time)
    plot!(eps, PCG_time)
    plot!(eps, ICCG_time)
    display(timings)

    return
end

# parameters for test
N = 64
h = 1/N
ϵ = 1e-6
opts = Dict("N"=>N,
            "h"=>h,
            "ϵ"=>ϵ,
            "rho"=>1.0
            )
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
# a bunch of test u
#actual_u = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
actual_u = [sin(8*x)+4*cos(8*y) for y in ys, x in xs]

# call tests
#visual_test(actual_u, xs, ys, opts)
timing_test(actual_u, opts)

