
using Makie, GLMakie
using Colors
using AbstractPlotting: textslider, colorswatch, Node, hbox, vbox
using AbstractPlotting

Makie.AbstractPlotting.inline!(false)


function showBoundaryCuboids_old(world::World)
    x = 1:world.dims[1]
    y = 1:world.dims[2]
    grid = zeros(world.dims[1], world.dims[2])
    for (idx_tup, value) in world.boundary_adjacency
        grid[idx_tup[1], idx_tup[2]] = 1
    end
    boundary_scene = meshscatter(
        x, y, grid,
    )
    parent = Scene(resolution = (1000, 800))
    interact = vbox(
        boundary_scene, parent = parent
    )
end

function showBoundaryCuboids(world::World)
    xs = 1:world.dims[1]
    ys = 1:world.dims[2]
    elev = world.dims[1]
    zs = 1:elev
    # needs to be a 1D array
    pts = [Point3f0(x, y, z) for x in xs for y in ys for z in zs]
    pointers = [Point3f0(0, 0, 0) for x in xs for y in ys for z in zs]
    # UPDATE here
    #   have cuboid at bidx point to aidx cuboid 
    #for (bidx, aidx) in world.boundary_adjacency
    #    change this indexing to index linearly, multiply, etc.
    #    pointers[bidx[1], bidx[2], elev] = Point1
    #end
    boundary_scene = arrows(
        pts, pointers 
    )
    parent = Scene(resolution = (1000, 800))
    interact = vbox(
        boundary_scene, parent = parent
    )
end