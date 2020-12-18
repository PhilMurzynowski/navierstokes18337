using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap
using BenchmarkTools

include("Matrices.jl")
include("CG.jl")

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
    hm, cl = get_u_scene(xs, ys, u)

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

    x_guess = zeros(N, N)
    tmp1 = zeros(N-2, N-2)
    tmp2 = copy(tmp1)
    tmp3 = copy(tmp1)
    tmp4 = copy(tmp1)
    # note, meant to use CG_Poisson here, but broken for some reason and PCG_Poisson_diag is equivalent
    CGPoisson_sol, iter = PCG_Poisson_diag(reshape(source, N, N), x_guess, ϵ, opts, tmp1, tmp2, tmp3, tmp4)
    hm_CGPoisson, cl_CGPoisson = get_u_scene(xs, ys, CGPoisson_sol)

    x_guess = zeros(N, N)
    tmp = copy(x_guess)
    Grid_sol, iter = grid_poisson_solve(x_guess, tmp, reshape(source, N, N)[2:end-1, 2:end-1], opts, 1e3)
    hm_grid, cl_grid = get_u_scene(xs, ys, Grid_sol)

    # hacky way to arrange things
    parent = Scene(resolution= (400, 800))
    full_scene = hbox(
                    vbox(hm_actual, cl_actual),
                    vbox(hm_CG, cl_CG), 
                    vbox(hm_PCG, cl_PCG), 
                    vbox(hm_ICCG, cl_ICCG), 
                    vbox(hm_CGPoisson, cl_CGPoisson), 
                    parent=parent)
    # separate plot for grid because its buggy
    parent = Scene(resolution= (300, 300))
    grid_scene = hbox(
                    vbox(hm_grid, cl_grid), 
                    parent=parent)
    #println(iter)
    #display(u)
    #display(Grid_sol)
    display(full_scene)
    display(grid_scene)
    return
end

"""
Given a source vector compare timings of different solution methods.
"""
function timing_test(source)
    return
end

#Test
#=
CG(A, b, x_guess, ϵ, max_iter=1e3)
PCG(A, Minv, b, x_guess, ϵ, max_iter=1e3)
CG_Poisson(b, x_guess, ϵ, opts, tmp1, tmp2, tmp3)
ICCG(A, U, b, x_guess, ϵ, max_iter=1e3)

Then retest all with BandedMatrices

# also test this! might need some reformatting
function grid_poisson_solve(p1, p2, velocity_component, opts, max_iter=500)
=#

# parameters for test
N = 32
h = 1/N
ϵ = 1e-4
opts = Dict("N"=>N,
            "h"=>h,
            "ϵ"=>ϵ,
            "rho"=>1.0
            )
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
actual_u = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
visual_test(actual_u, xs, ys, opts)
#=
actual_pressure = zeros(N, N)
#actual_pressure = [x*y for y in ys, x in xs]
#actual_pressure = [sin(2*x)+4*cos(4*y) for y in ys, x in xs]
=#

#hm_actual = heatmap(xs_extend, ys_extend, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
#cl_actual = colorlegend(hm_actual[end], raw = true, camera = campixel!)

