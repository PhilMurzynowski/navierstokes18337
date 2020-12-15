using LinearAlgebra, Kronecker
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap

include("CG.jl")

# methods dedicated to solving Poisson Eq

function genPoissonMtx(opts)

    N = opts["N"]
    hi = 1/opts["h"]
    println("updated gen mtx")

    L = zeros(N, N)
    L[diagind(L, 0)] .= 2*hi^2
    L[diagind(L, -1)] .= -1*hi^2
    L[diagind(L, 1)] .= -1*hi^2
    # build Poisson pressure matrix with Kronecker
    # cant use I with Kronecker, not working
    Id = 1*Matrix(I, N, N)
    P = (Id ⊗ L) + (L ⊗ Id)
    # Von Neumann BC and reference point modifications
    BC = zeros(N, N)
    BC[1, 1] = -hi^2
    BC[N, N] = -hi^2
    P += (Id ⊗ BC)
    P[1:N, 1:N] += -hi^2*Id
    P[N*(N-1)+1:N^2, N*(N-1)+1:N^2] += -hi^2*Id
    #P[1, :] .= 0
    #P[1, 1] = 1*hi^2
    #P[N^2, :] .= 0
    #P[N^2, 1] = 1*hi^2
    return P
end

# test setup
# generate pressure field and determine source terms
# then try to go back backwards by solving Poissson Eq with source terms

N = 64
h = 0.01
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
#actual_pressure = [x*y for y in ys, x in xs]
#actual_pressure = [sin(2*x)+4*cos(4*y) for y in ys, x in xs]
actual_pressure = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
#display(actual_pressure)
# not showing axis because top left corner of matrix is meant to be plotted in top left
# but then axis put bottom left as origin (understandably, which makes things confusing)

# can make these surface plots instead later! would be cooler!
# heatmap also low res at the moment

hm_actual = heatmap(xs, ys, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
title("Actual Pressure")
cl_actual = colorlegend(hm_actual[end], raw = true, camera = campixel!)



opts = Dict("N"=>N,
            "Nx"=>N,
            "Ny"=>N,
            "dx"=>h,
            "dy"=>h,
            "h"=>h,
            )

P = genPoissonMtx(opts)
b = P*vec(actual_pressure)
x_guess = zeros(length(b))
ϵ = 1e-10
pressure_CG_std, num_iter_CG_std = CG_std(P, b, x_guess, ϵ)
#display(pressure_CG_std)
pressure_CG_std = reshape(pressure_CG_std, N, N)

hm_CG_std = heatmap(xs, ys, reverse(pressure_CG_std', dims=2), colormap=:berlin, show_axis=false)
cl_CG_std = colorlegend(hm_CG_std[end], raw = true, camera = campixel!)
#title(full_scene, "hi")
#display(full_scene)


println("hi")
pressure_CG_Poisson, num_iter_CG_Poisson = CG_Poisson(b, zeros(N+2, N+2), ϵ, opts)
xs_extend = 0.0:h:h*(N+1)
ys_extend = 0.0:h:h*(N+1)
hm_CG_Poisson = heatmap(xs_extend, ys_extend, reverse(pressure_CG_Poisson', dims=2), colormap=:berlin, show_axis=false)
cl_CG_Poisson = colorlegend(hm_CG_Poisson[end], raw = true, camera = campixel!)

# plotting
# hacky way to arrange things
# can arrange things more cleverly if dont look good
parent = Scene(resolution= (1000, 500))
full_scene = hbox(
                vbox(hm_actual, cl_actual),
                vbox(hm_CG_std, cl_CG_std), 
                vbox(hm_CG_Poisson, cl_CG_Poisson), 
                parent=parent)
println("bye")
display(full_scene)
display(P)