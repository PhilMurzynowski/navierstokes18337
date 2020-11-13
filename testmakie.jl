using Makie, AbstractPlotting


# going through Makie tutorial

Makie.AbstractPlotting.inline!(true)
scene = Scene()
display(scene)
#points = [Point{2, Float32}((cos(t), sin(t))) for t in LinRange(0, 2pi, 20)]
#colors = UnitRange{Int32}(1:20)
#scatter!(scene, points, color = colors, markersize = 10, show_axis = false)
points = [Point2f0((cos(t), sin(t))) for t in LinRange(0, 2pi, 20)]
colors = 1:20
scatter!(scene, points, color = colors, markersize = 10)
scatterobject = scene[end]
scatterobject.markersize = 30
display(scene)

# using Observables
xs = -pi:0.01:pi
frequency = Node(3.0)   # Node is Observable
phase = Node(0.0)
# lift creates new observable by applying function to observables
# do creates an anonymous function and passes it as argument
ys = lift(frequency, phase) do fr, ph
    # believe @. macro say to perform elementwise operation
    @. 0.3 * sin(fr * xs - ph)
end

# :blue is a symbol named blue
lines!(scene, xs, ys, color = :blue, linewidth = 3)

# observables mutated with empty brackets
frequency[] = 20
# now hook up with sliders and buttons!
save("project/tutorial_sineplot.png", scene)

# creating animations
framerate = 30
timestamps = 0:1/framerate:3
record(scene, "project/tutorial_phase_animation.mp4", timestamps, framerate = framerate) do t
    phase[] = 2 * t * 2pi
end

# only one sine wave moving as accidentally called lines! again, did not extract object from scene

# 3D plotting
scene3D = Scene(camera = cam3d!)
cam = cameracontrols(scene3D)
cam.rotationspeed[] = 0.05
helix = [Point3f0((cos(10*t), sin(10*t), t)) for t in LinRange(0, 2pi, 1000)]
scatter!(scene3D, helix, color = colors, markersize = 50)

# trying to get interactive 3D
# copied from https://nextjournal.com/julia-berlin/makie-demo 
# removed all Interact as that seems to only work for Web based demos
using Colors, AbstractPlotting, GLMakie
AbstractPlotting.set_theme!() # reset theme
function lorenz(t0, a, b, c, h)
     Point3f0(
         t0[1] + h * a * (t0[2] - t0[1]),
         t0[2] + h * (t0[1] * (b - t0[3]) - t0[2]),
         t0[3] + h * (t0[1] * t0[2] - c * t0[3]),
     )
end
 # step through the `time`
 function lorenz(array::Vector, a = 5.0 ,b = 2.0, c = 6.0, d = 0.01)
     t0 = Point3f0(0.1, 0, 0)
     for i = eachindex(array)
         t0 = lorenz(t0, a,b,c,d)
         array[i] = t0
     end
     array
 end
mcolor = widget(colorant"red")
n1, n2 = 18, 30
N = n1*n2
args_n = observe.((a, b, c, d))
v0 = lorenz(zeros(Point3f0, N), to_value.(args_n)...)
positions = lift(lorenz, Makie.Node(v0), args_n...)
mscene = meshscatter(
  positions,
  markersize = scales,
  intensity = collect(range(0f0, stop = 1f0, length = length(positions[]))),
  color = mcolor, center = false
)

# this setup below allows for interactive rotation
# https://makie.juliaplots.org/dev/plotting_functions.html#Plotting-Functions
using GLMakie
using AbstractPlotting
xs = cos.(1:0.5:20)
ys = sin.(1:0.5:20)
zs = LinRange(0, 3, length(xs))
meshscatter(xs, ys, zs, markersize = 0.1, color = zs)

# I am likely more interested in volume, heatmap or just scatter
using GLMakie
using AbstractPlotting
r = LinRange(-1, 1, 100)
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
cube_with_holes = cube .* (cube .> 1.4)
volume(cube_with_holes, algorithm = :iso, isorange = 0.05, isovalue = 1.7)



#=

Trying to add sliders 

If you need a normal Makie scene in a layout,
for example for 3D plots, you have to use LScene right now.
It's just a wrapper around the normal Scene that makes it layoutable.
The underlying Scene is accessible via the scene field.
You can plot into the LScene directly, though.

You can pass keyword arguments to the underlying Scene object
to the scenekw keyword. Currently, it can be necessary to pass a
couple of attributes explicitly to make sure they are not inherited
from the main scene (which has a pixel camera and no axis, e.g.).

=#
# to see this visualization just type scene into REPL or display(scene)
# works, no sliders yet, but LScene functional
using Makie
using AbstractPlotting.MakieLayout
scene, layout = layoutscene(resolution = (1400, 900))
lscenes = layout[1:2, 1:3] = [LScene(scene, camera = cam3d!, raw = false) for _ in 1:6]
[scatter!(lscenes[i], rand(100, 3), color = c)
    for (i, c) in enumerate([:red, :blue, :green, :orange, :black, :gray])]
sl1 = LSlider(scene, range = 0:0.01:10, startvalue = 3, halign = :right, valign = :bottom)
sl2 = LSlider(scene, range = 0:0.01:10, startvalue = 5, halign = :right, valign = :top)

# not working

using GLMakie
using AbstractPlotting
using AbstractPlotting.MakieLayout
scene, layout = layoutscene(resolution = (1200, 900))
lscene = layout[1, 1] = LScene(scene, scenekw = (camera = cam3d!, raw = true))
r = LinRange(-1, 1, 100)
#sl1 = LSlider(scene, range = 0:0.01:10, star
#sl2 = LSlider(scene, range = 0:0.01:10, startvalue = 5, halign = :right, valign = :top)
cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
cube_with_holes = cube .* (cube .> 1.4)
cube_with_holes = cube .* (cube .> cut_threshold)
volume!(lscene, cube_with_holes, algorithm = :iso, isorange = 0.05, isovalue = 1.7)

# more examples for slider functionality
# 2D, not sure why have s2[end] and s1[end]
#   might be similars to observables needed to be followed by [] (indexing syntax)
# http://juliaplots.org/MakieReferenceImages/gallery//sliders/index.html
# 3D example!
# http://juliaplots.org/MakieReferenceImages/gallery//volume_slider/index.html
# but they seem outdated and need to be modified for new infrastructure

# start with 2D example, presumably easier to get running
# got it to plot in 3D actually!
# sliders are not fully functional, clipped, likely not entirely correct values, but still really good
using GLMakie, AbstractPlotting
scene, layout = layoutscene(resolution = (1200, 900))
lscene = layout[1, 1] = LScene(scene, scenekw = (camera = cam3d!, raw = true))
sl1 = LSlider(scene, range = 0:0.01:10, startvalue = 3, halign = :right, valign = :bottom)
sl2 = LSlider(scene, range = 40:1:100, startvalue = 40, halign = :right, valign = :top)
data = lift(sl1.value) do v
    map(LinRange(0, 2pi, 100)) do x
        4f0 .* Point2f0(sin(x) + (sin(x * v) .* 0.1), cos(x) + (cos(x * v) .* 0.1))
    end
end
scatter!(lscene, data, markersize = sl2.value)

# 3D example
# getting Realtime updates but glitches a lot
# disappers if zoom out too much, starts super zoomed in

using GLMakie, AbstractPlotting
scene, layout = layoutscene(resolution = (1200, 900))
lscene = layout[1, 1] = LScene(scene, scenekw = (camera = cam3d!, raw = true))
sl1 = LSlider(scene, range = 1:0.01:2, startvalue = 1.4, halign = :right, valign = :bottom)
r = LinRange(-1, 1, 100)
data = lift(sl1.value) do cut_threshold
    cube = [(x.^2 + y.^2 + z.^2) for x = r, y = r, z = r]
    cube_with_holes = cube .* (cube .> cut_threshold)
end
volume!(lscene, data, algorithm = :iso, isorange = 0.05, isovalue = 1.7)

# another attempt
# A perfect working example!!!!
Makie.AbstractPlotting.inline!(false)
using Colors
using AbstractPlotting: textslider, colorswatch, Node, hbox, vbox
using AbstractPlotting

# why does textslider return two things? not aliasing probably,
# something like an object and value, but more complicated with linking?

s1, a = textslider(0f0:50f0, "a", start = 13)
s2, b = textslider(-20f0:20f0, "b", start = 10)
s3, c = textslider(0f0:20f0, "c", start = 2)
s4, d = textslider(range(0.0, stop = 0.02, length = 100), "d", start = 0.01)
s5, scales = textslider(range(0.01, stop = 0.5, length = 100), "scale", start = 0.1)
s6, colorsw, pop = colorswatch()

function lorenz(t0, a, b, c, h)
    Point3f0(
        t0[1] + h * a * (t0[2] - t0[1]),
        t0[2] + h * (t0[1] * (b - t0[3]) - t0[2]),
        t0[3] + h * (t0[1] * t0[2] - c * t0[3]),
    )
end
# step through the `time`
function lorenz(array::Vector, a = 5.0 ,b = 2.0, c = 6.0, d = 0.01)
    t0 = Point3f0(0.1, 0, 0)
    for i = eachindex(array)
        t0 = lorenz(t0, a,b,c,d)
        array[i] = t0
    end
    array
end
n1, n2 = 18, 30
N = n1*n2
args_n = (a, b, c, d)
v0 = lorenz(zeros(Point3f0, N), to_value.(args_n)...)
positions = lift(lorenz, Node(v0), args_n...)
# not entirely clear on rotations
rotations = lift(diff, positions)
rotations = lift(x-> push!(x, x[end]), rotations)

mesh_scene = meshscatter(
    positions,
    markersize = scales,
    # what is rotation here?
    rotation = rotations,
    # not sure what intensity is either
    intensity = collect(range(0f0, stop = 1f0, length = length(positions[]))),
    color = colorsw
)
# not sure about what this parent, vbox, and hbox structures are either
parent = Scene(resolution = (1000, 800))
interact = vbox(
    hbox(s1, s2, s3, s4, s5, s6),
    mesh_scene, parent = parent
)