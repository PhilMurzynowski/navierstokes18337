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

const BoundaryAdjacency = Dict{Tuple{Int, Int}, Tuple{Int, Int}}

"""
This seems to OOP, refactor later.
Currently assuming worlds also have cuboid shape with no obstacles within.
Cuboids on periphery are boundary cuboids, are not updated in the regular fashion.
    They model boundary conditions.
"""
mutable struct World{T, N}
    dims::Array{Int, 1}
    boundary_adjacency::BoundaryAdjacency
    cuboids::Array{Cuboid{T}, N}
    
    # Intialize all cells 
    # dimensions is size in y, x, respectively
    function World(dimension::Int, dimensions::Array{Int, 1}, type=Float64)
        if dimension == 2
            cuboids = Array{Cuboid{type}, dimension}(undef, dimensions[1], dimensions[2])
            cuboid = Cuboid(zeros(2), 0.0)
            # initialize boundary adjacency for all blocks on perimiter
            # as assuming cuboid world with no internal obstacles
            boundary_adjacency = BoundaryAdjacency()
            # corners
            push!(boundary_adjacency, (1, 1) => (2, 2))
            push!(boundary_adjacency, (dimensions[1], 1) => (dimensions[1]-1, 2))
            push!(boundary_adjacency, (1, dimensions[2]) => (2, dimensions[2]-1))
            push!(boundary_adjacency, (dimensions[1], dimensions[2]) => (dimensions[1]-1, dimensions[2]-1))
            for j in 2:dimensions[1]-1 
                # left col
                push!(boundary_adjacency, (j, 1) => (j, 2))
                # right col
                push!(boundary_adjacency, (j, dimensions[2]) => (j, dimensions[2]-1))
            end
            for i in 2:dimensions[2]-1
                # top row
                push!(boundary_adjacency, (1, i) => (2, i))
                # bottom row
                push!(boundary_adjacency, (dimensions[1], i) => (dimensions[1]-1, i))
            end
        elseif dimension == 3    
            cuboids = Array{Cuboid{type}, dimension}(undef, dimensions[1], dimensions[2], dimensions[3])
            cuboid = Cuboid(zeros(3), 0.0)
        end
        fill!(cuboids, cuboid)
        new{type, dimension}(dimensions, boundary_adjacency, cuboids)
    end
end
