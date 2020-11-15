module pg
include("Structures.jl")
include("Solvers.jl")
include("World.jl")
include("Plotting.jl")
world = World(2, [100, 100])
interact = showBoundaryCuboids(world)
end