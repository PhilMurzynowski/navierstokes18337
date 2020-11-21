module pg
include("Structures.jl")
include("Solvers.jl")
include("World.jl")
include("Plotting.jl")
world = World(1.0, 1.0, 1.0, 2, [10, 10])
world2 = deepcopy(world)
interact = showBoundaryCuboids(world)
# set lid velocity
# to both worlds
for i in 1:world.dims[2]-1
    world.cuboids[1, i].u[1] = 1.0
end
for i in 1:world.dims[2]-1
    world2.cuboids[1, i].u[1] = 1.0
end
simple_nomtx(world, world2)
end