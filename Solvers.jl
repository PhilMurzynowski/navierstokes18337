#=
Solver
=#

# does not make use of matrix operations
# assuming world is 2D
# https://arxiv.org/ftp/physics/papers/0407/0407002.pdf
# ignore this one for now https://www.montana.edu/mowkes/research/source-codes/GuideToCFD.pdf 
function simple_nomtx(world:World)
    dt, dx, dy = world.deltas
    # naive iteration, not optimized for min cache misses
    # add in check for boundary conditions
    for j in 1:world.dims[1]
        for i in 1:world.dims[2]
            # each cell has left and bottom
            # j is increasing going down
            # i is increasing going right

            # really messy naming, cleanup
            u_left = world.cells[j, i].u[1]
            v_bottom = world.cells[j, i].u[2]
            u_right = world.cells[j, i+1].u[1]
            v_top = world.cells[j-1, i].u[2]
            u_next_right = world.cells[j, i+2].u[1]
            v_next_bottom = world.cells[j+1, i].u[2]
            u_right_down = world.cells[j+1, i+1].u[1] 
            u_right_up = world.cells[j-1, i+1].u[1] 
            v_bottom_right = world.cells[j, i+1].u[2] 
            v_bottom_left = world.cells[j, i-1].u[2] 
            v_top_right = world.cells[j-1, i+1].u[2]
            v_right = 1/4*(v_bottom + v_bottom_right + v_top + v_top_right)

            a3 = (u_next_right - 2*u_right + u_left) / dx^2
            a4 = (u_right_down - 2*u_right + u_right_up) / dy^2
            b3 = (v_next_bottom - 2*v_bottom + v_top) / dy^2
            b4 = (v_bottom_right - 2*v_bottom + v_bottom_left) / dx^2
            a1 = -(u_right^2 - u_left^2) / (2*dx) - () / (2*dy)
        end
    end
end

function piso(world::World)
    return
end

# http://webx.ubi.pt/~pjpo/ri23.pdf
function piso_improved(world::World)
    return
end
