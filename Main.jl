using Profile

"""
Main file for entire directory.
Also serves as example for how to run a simulation.

1. Initialize a problem
2. Pass flags for setup and solving scheme
"""

include("gridSolver.jl")
include("VorticityStream.jl")

GLMakie.activate!()
Makie.AbstractPlotting.inline!(false)

# define velocity boundary conditions at the four walls
BC_opts = Dict("u_top"=>1.0,
               "u_bottom"=>0.0,
               "v_left"=>0.0,
               "v_right"=>0.0)

# define general options and specific parameters
# e.g. General: timestep size, simulation size
# but for VorticityStream Solver may want flags for how to solve Poisson Eq
# e.g whether to us CG, PCG, ICCG, specified in CG.jl

N = 32
h = 1/N         # chose this spacing largely for visuals
dt = 0.001      # timestep
opts = Dict("timesteps"=>400,   # number of steps to simulate                
            "ϵ"=>1e-8,          # tolerance for Solving Poisson Eqauation       
            "N"=>N,
            "Nx"=>N,            # Nx, Ny left in here equal to N as a reminder that the region is square
            "Ny"=>N,            # admittedly clunky, left in for update to non-square regions
            "dx"=>h,            # must have dx=dy as well
            "dy"=>h,
            "h"=>h,
            "dt"=>dt,
            "rho"=>1.0,         # ρ only available in GridSolver currently, so keep at 1 for comparison
            "Re"=>10,           # Reynolds number, easier to see simulation with lower Reynolds number
            )

# run simulation with desired solver
# choice 1: Grid
# choice 2: VorticityStream
#solver = "Grid"
#solver = "VorticityStream"
solver = "VorticityStream"

if solver == "Grid"
    u, v, p = run_grid_simulation(opts, BC_opts)
    plot_uvp(u, v, p, opts)
elseif solver == "VorticityStream"
    ω, ψ = run_vorticitystream_simulation(opts, BC_opts)
    plot_ωψ(ω, ψ)
else
    # runs both without plotting
    @time u, v, p = run_grid_simulation(opts, BC_opts)
    @time ω, ψ = run_vorticitystream_simulation(opts, BC_opts)
    #@profile u, v, p = run_grid_simulation(opts, BC_opts)
    #Profile.print()
    #@profile ω, ψ = run_vorticitystream_simulation(opts, BC_opts)
    #Profile.print()
    return
end