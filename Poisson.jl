using Makie, AbstractPlotting
using AbstractPlotting: Node, hbox, vbox, heatmap

# methods dedicated to solving Poisson Eq

# test setup
# generate pressure field and determine source terms
# then try to go back backwards by solving Poissson Eq with source terms
#correct_pressure = [sin(x)+cos(y) for x in 0.0:1:1024 for y in 0.0:1:1024]

N = 1024
h = 0.01
xs = 0.0:h:h*(N-1)
ys = 0.0:h:h*(N-1)
#actual_pressure = [x*y for y in ys, x in xs]
actual_pressure = [sin(2*x)+4*cos(4*y) for y in ys, x in xs]
#actual_pressure = [sinc(sqrt(x^2+y^2)) for y in ys, x in xs]
display(actual_pressure)
# not showing axis because top left corner of matrix is meant to be plotted in top left
# but then axis put bottom left as origin (understandably, which makes things confusing)

# can make these surface plots instead later! would be cooler!

hm = heatmap(xs, ys, reverse(actual_pressure', dims=2), colormap=:berlin, show_axis=false)
cl = colorlegend(hm[end], raw = true, camera = campixel!)
parent = Scene(resolution= (1000, 500))
full_scene = vbox(vbox(hm, cl), parent=parent)
display(full_scene)