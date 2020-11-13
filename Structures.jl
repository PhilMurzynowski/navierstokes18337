#=
Grid discretization primitives.
=#

# considering using abstract types
#abstract type Cuboid end

#=
Staggered grid
    u[1] is wrt west face
    u[2] is wrt south face
    u[3] is wrt to bottom face
    pressure is wrt center
=#
mutable struct Cuboid{T}
    u::Array{T, 1}
    pressure::T
end


mutable struct World{T, N}
    cells::Array{Cuboid{T}, N}
    
    # Intialize all cells 
    function World(dimension, dimensions, type=Float64)
        if dimension == 2
            world = Array{Cuboid{type}, dimension}(undef, dimensions[1], dimensions[2])
            cuboid = Cuboid(zeros(2), 0.0)
        elseif dimension == 3    
            world = Array{Cuboid{type}, dimension}(undef, dimensions[1], dimensions[2], dimensions[3])
            cuboid = Cuboid(zeros(3), 0.0)
        end
        fill!(world, cuboid)
        new{type, dimension}(world)
    end
end