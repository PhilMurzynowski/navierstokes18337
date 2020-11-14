#=
Grid discretization primitives.
=#

# considering using abstract types
#abstract type Cuboid end

#=
Staggered grid for velocity
Pressure is wrt center
=#
mutable struct Cuboid{T}
    u::Array{T, 1}
    pressure::T
end


mutable struct World{T, N}
    dims::Array{Float32, 1}
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
        new{type, dimension}(dimensions, world)
    end
end