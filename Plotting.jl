
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

# function to take in matrices of U and V and output streamplots
#function streamplotUV(u, v)
#    idx_into_uv(x::Point2{T}, t) where T = Point2{T}(u[x[1]], v[x[2]])
#
#    sf = Node(Base.Fix2(idx_into_uv, 0e0))
#    
#    title_str = Node("t = 0.00")
#    
#    sp = streamplot(
#            sf,
#            -2..2, -2..2;
#            linewidth = 2,
#            padding = (0, 0),
#            arrow_size = 0.09,
#            colormap =:magma
#        )
#    
#    sc = title(sp, title_str)
#    
#    record(sc, "output.mp4", LinRange(0, 20, 5*30)) do i
#      sf[] = Base.Fix2(v, i)
#      title_str[] = "t = $(round(i; sigdigits = 2))"
#    end
#end