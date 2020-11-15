module pg
include("Structures.jl")
include("Solvers.jl")
include("World.jl")
include("Plotting.jl")
world = World(2, [20, 20])
interact = showBoundaryCuboids(world)
end