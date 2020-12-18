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

N = 3
h = 0.01
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
xs_extend = 0.0:h:h*(N+1)
ys_extend = 0.0:h:h*(N+1)
#actual_pressure = [x*y for y in ys, x in xs]
#actual_pressure = [sin(2*x)+4*cos(4*y) for y in ys, x in xs]
pressure = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
# make it N+2 by N+2 for BC
actual_pressure = zeros(N+2, N+2)
actual_pressure[2:end-1, 2:end-1] = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
actual_pressure[:, 1] = actual_pressure[:, 2]
actual_pressure[:, end] = actual_pressure[:, end-1]
actual_pressure[1, :] = actual_pressure[2, :]
actual_pressure[end, :] = actual_pressure[end-1, :]
# should it be 0?
#actual_pressure[end, :] .= 0                   # do this last to avoid overwriting
#display(actual_pressure)
# not showing axis because top left corner of matrix is meant to be plotted in top left
# but then axis put bottom left as origin (understandably, which makes things confusing)

# can make these surface plots instead later! would be cooler!
# heatmap also low res at the moment

#hm_actual = heatmap(xs, ys, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
hm_actual = heatmap(xs_extend, ys_extend, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
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
b = P*vec(@view actual_pressure[2:end-1, 2:end-1])

#=
x_guess = zeros(length(b))
ϵ = 1e-10
# can't really handle BC that well it seems
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

=#

function PoissonMult(x, opts)
    N = opts["N"]
    Δx = opts["dx"]
    Δy = opts["dy"]
    # reshaped for readability
    #input = reshape(input, length(input), 1)
    input = reshape(input, N, N)

    # handle top left square edge case
    # first N entries

    # outermost loop, loop over N rows in Laplacian at a time
    for k in 2:N-1
        for i in 2:N+1
            for j in 2:N+1
                output[j, i] = 1/Δx^2*(input[j, i+1] - 2*input[j, i] + input[j, i-1])
                output[j, i] += 1/Δy^2*(input[j+1, i] - 2*input[j, i] + input[j-1, i])
            end
        end
    end

    # handle bottom right square edge case
    # last N entries
end

#= scrap
        for i in 2:N+1
            for j in 2:N+1
                output[j, i] = 1/Δx^2*(input[j, i+1] - 2*input[j, i] + input[j, i-1])
                output[j, i] += 1/Δy^2*(input[j+1, i] - 2*input[j, i] + input[j-1, i])
            end
        end
=#

# test Poisson Mtx mult
println("testing mult")
display(P)
println("normal mult")
display(P*vec(pressure))
println("special mult")