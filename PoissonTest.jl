using LinearAlgebra
using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap

include("Matrices.jl")
include("CG.jl")

"""
File dedicated to testing and benchmarking different approaches to solving Poisson Eq.
Testing approach:
    1. A pressure field is generated
    2. Source terms are determined by multiplying pressure field by laplacian
    3. Validate solver by checking pressure field after solving Poisson Eq given source terms
Note: A nice update would be to use surface plots instead later! would be cooler!
"""

"""
Given a pressure field to test against, return the appropriate source vector
"""
function determine_source(pressure_field, size, spacing)
    P = genPoissonMtx(size, spacing)
    source = P*vec(pressure_field)
    return source
end

#=
opts = Dict("N"=>N,
            "h"=>h,
            )
=#
# parameters for test
N = 32
h = 1/N
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
xs_extend = 0.0:h:h*(N+1)
ys_extend = 0.0:h:h*(N+1)
actual_pressure = zeros(N, N)
#actual_pressure = [x*y for y in ys, x in xs]
#actual_pressure = [sin(2*x)+4*cos(4*y) for y in ys, x in xs]
actual_pressure = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]

hm_actual = heatmap(xs_extend, ys_extend, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
cl_actual = colorlegend(hm_actual[end], raw = true, camera = campixel!)

# plotting
# hacky way to arrange things
# can arrange things more cleverly if dont look good
parent = Scene(resolution= (1000, 500))
full_scene = hbox(
                vbox(hm_actual, cl_actual),
                #vbox(hm_CG_std, cl_CG_std), 
                #vbox(hm_CG_Poisson, cl_CG_Poisson), 
                parent=parent)
display(full_scene)
