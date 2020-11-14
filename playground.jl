module pg

include("Structures.jl")
include("Piso.jl")
include("World.jl")
world = World(2, [10, 10])

end