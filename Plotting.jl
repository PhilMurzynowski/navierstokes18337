
using Makie, GLMakie
using Colors
using AbstractPlotting: textslider, colorswatch, Node, hbox, vbox
using AbstractPlotting

Makie.AbstractPlotting.inline!(false)

function showBoundaryCuboids(world::World)
    xmax = world.dims[1]
    ymax = world.dims[2]
    elev = 5
    zmax = elev
    # needs to be a 1D array
    pts = [Point3f0(x, y, z) for x in 1:xmax for y in 1:ymax for z in 1:zmax]
    pointers = [Point3f0(0, 0, 0) for x in 1:xmax for y in 1:ymax for z in 1:zmax]
    for (bidx, aidx) in world.boundary_adjacency
        pointers[(bidx[1]-1)*ymax*zmax + (bidx[2]-1)*zmax + elev] = Point3f0(aidx[1]-bidx[1], aidx[2]-bidx[2], -elev)
    end
    boundary_scene = arrows(
        pts, pointers 
    )
    parent = Scene(resolution = (1000, 800))
    interact = vbox(
        boundary_scene, parent = parent
    )
end