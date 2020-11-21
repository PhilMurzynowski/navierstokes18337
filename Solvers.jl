#=
Solver
=#

# does not make use of matrix operations
# assuming world is 2D
# https://arxiv.org/ftp/physics/papers/0407/0407002.pdf
# ignore this one for now https://www.montana.edu/mowkes/research/source-codes/GuideToCFD.pdf 

"""
Function to check whether indeces for a given world refer to a 
cuboid on a boundary.
Intended for debugging.
"""
function boundaryCuboid(world::World, row::Int, col::Int)
    if (row == 1) || (col == 1) || (row == world.dims[1]) || (col == world.dims[2])
        return true        
    end
    return false
end

#=
"""
Update cuboid velocity at index [j, i].
Likely will not be used as requires branch for every index,
which is unnecessary.
"""
"""
function updateCuboidVelocity(world::World, new_world::World, j::Int, i::Int)
    if boundaryCuboid(world, j, i)
        updateBoundaryCuboidVelocity(world, new_world, j, i)
    else
        updateCenterCuboidVelocity(world, new_world, j, i)
    end
end
"""
=#

"""
Update the velocity of the cuboid at index [j, i].
The cuboid should not be a boundary cuboid of the world.
"""
function predictCuboidVelocity(world::World, new_world::World, j::Int, i::Int)
    Δt = world.Δt
    Δx, Δy = world.Δxyz, world.Δxyz
    # 1 / Re
    Re_inv = world.Re_inv
    # each cell has right and bottom
    # j is increasing going down
    # i is increasing going right
    # really messy naming, cleanup
    # u's
    u_left = world.cuboids[j, i-1].u[1]
    u_right = world.cuboids[j, i].u[1]
    u_twice_right = world.cuboids[j, i+1].u[1]
    u_right_down = world.cuboids[j+1, i].u[1] 
    u_left_down = world.cuboids[j+1, i-1].u[1] 
    u_right_up = world.cuboids[j-1, i].u[1] 
    # v's
    v_bottom = world.cuboids[j, i].u[2]
    v_top = world.cuboids[j-1, i].u[2]
    v_twice_bottom = world.cuboids[j+1, i].u[2]
    v_bottom_right = world.cuboids[j, i+1].u[2] 
    v_bottom_left = world.cuboids[j, i-1].u[2] 
    v_top_right = world.cuboids[j-1, i+1].u[2]
    # linearly interpolate (average)
    u_bottom = 1/4*(u_right + u_left + u_right_down + u_left_down)
    v_right = 1/4*(v_bottom + v_bottom_right + v_top + v_top_right)
    # group together in stages for readability
    a3 = (u_twice_right - 2*u_right + u_left) / Δx^2
    a4 = (u_right_down - 2*u_right + u_right_up) / Δy^2
    b3 = (v_twice_bottom - 2*v_bottom + v_top) / Δy^2
    b4 = (v_bottom_right - 2*v_bottom + v_bottom_left) / Δx^2
    # calculations using interpolated velocites
    a1 = -(u_twice_right^2 - u_left^2)/(2*Δx) - v_right*(u_right_down - u_right_up)/(2*Δy)
    b1 = -(v_twice_bottom^2 - v_top^2)/(2*Δy) - u_bottom*(v_bottom_right - v_bottom_left)/(2*Δx)
    # believe paper wrongly negated a1 here, check
    A = a1 + Re_inv*(a3 + a4)
    B = b1 + Re_inv*(b3 + b4)
    # u^{n+1}
    # LHS of eqns 7 & 8
    p_center = world.cuboids[j, i].pressure
    p_right = world.cuboids[j, i+1].pressure
    p_down = world.cuboids[j+1, i].pressure
    u_right_next = u_right + Δt*(A - (p_right - p_center)/Δx)
    v_bottom_next = v_bottom + Δt*(B - (p_down - p_center)/Δy)
    # update
    new_world.cuboids[j, i].u = [u_right_next, v_bottom_next]
end

"""
Simple method to solve for pressure correction
"""
function calcPressureCorrection(world::World, new_world::World, j::Int, i::Int)
    Δt = world.Δt
    Δx, Δy = world.Δxyz, world.Δxyz
    p_right = world.cuboids[j, i+1].pressure
    p_left = world.cuboids[j, i-1].pressure
    p_up = world.cuboids[j-1, i].pressure
    p_down = world.cuboids[j+1, i].pressure
    u_right = world.cuboids[j, i].u[1]
    u_left = world.cuboids[j, i-1].u[1]
    v_bottom = world.cuboids[j, i].u[2]
    v_top = world.cuboids[j-1, i].u[2]
    a = 2*Δt*(1/Δx^2 + 1/Δy^2)
    b = -Δt/Δx^2
    c = -Δt/Δy^2
    d = 1/Δx*(u_right - u_left) + 1/Δy*(v_bottom - v_top)
    return -1/a*(b*(p_right + p_left) + c*(p_down + p_up) + d)
end

"""
Simple method for velocity correction
"""
function calcVelocityCorrection(world::World, Δp, Δp_right, Δp_down)
    Δt = world.Δt
    Δx, Δy = world.Δxyz, world.Δxyz
    Δu = 1 / (Δt*Δx) * (Δp_right - Δp)
    Δv = 1 / (Δt*Δy) * (Δp_down - Δp)
    return Δu, Δv
end

"""
If want to implement a boundary condition such as dp/dx (pressure differential)
at a wall is 0 then need to update the boundary cubes with the pressure of the cube
nearest to it.

Keep hashmap of closest cube to each boundary cube for this purpose.
Hashamp is world.boundary_adjacency
    bidx => aidx
Where bidx is the key for corresponding to the indeces for each boundary cuboid,
and aidx is the value corresponding to the indeces for the cuboid marked as adjacent
to the boundary cuboid.

Simply iterating over all boundary cuboids should result in large upper bound on cache misses.
E.g. leftmost and rightmost boundary cuboids will have a cache miss for every access for large n.
"""
function updateBoundaryCuboidsPressure(world::World)
    for (bidx, aidx) in world.boundary_adjacency
       world.cuboids[bidx[1], bidx[2]].pressure = world.cuboids[aidx[1], aidx[2]].pressure
    end
end

"""
Update velocity values of entire grid.
Naive iteration, not optimized for locality and cache misses.
"""
function simple_nomtx(world::World, new_world::World)
    # assume boundary cuboids on periphery
    # loop over interior cuboids
    for j in 2:world.dims[1]-1
        for i in 2:world.dims[2]-1
            predictCuboidVelocity(world, new_world, j, i)
        end
    end
    for j in 2:world.dims[1]-1
        for i in 2:world.dims[2]-1
            calcPressureCorrection(world, new_world, j, i)
        end
    end
    # should also get dp for boundary cuboids

end

function piso(world::World)
    return
end

# http://webx.ubi.pt/~pjpo/ri23.pdf
function piso_improved(world::World)
    return
end
