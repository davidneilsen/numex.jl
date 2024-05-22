#=-----------------------------------------------------------------------
 =
 = Wed May 22 08:58:58 MDT 2024
 =
 = This file defines the primary data structures for the code.  
 =
 = The Grid data structure stores the spatial grid and information
 = about the grid.  The GridFunction data structure stores the functions
 = that are defined on the grid.
 =
 =-----------------------------------------------------------------------=#


#=-----------------------------------------------------------------------
 = 
 =    The grid points are in the array x.
 =
 =    nx = the number of points in x
 =    dx = the uniform spacing between points, dx = x_{i+1} - x_{i} for any i
 =    [xmin, xmax] = the range of the grid
 =    dt = the time step size, computed using the CFL condition, dt = CFL * dx.
 = 
 =-----------------------------------------------------------------------=#
struct Grid

    xmin :: Float64
    xmax :: Float64
    nx :: Int64
    x :: Vector{Float64}
    dx :: Float64
    dt :: Float64

    function Grid( nx0, bbox, cfl)

        nx = nx0
        xmin = bbox[1]
        xmax = bbox[2]
        dx = (xmax - xmin)/ (nx - 1)
        c = 1.0
        dt = cfl * dx
        x = LinRange(xmin, xmax, nx)
        new( xmin, xmax, nx, x, dx, dt)

    end

end

#=-----------------------------------------------------------------------
 = 
 =    All functions defined on the grid.
 =
 =    grid = the grid structure, defined above
 =    neqs = the number of evolved variables in the set of PDEs
 =    u = a state vector of evolved variables at the current time.
 =    u2 = work space for intermediate values of u used by an RK routine.
 =    dxu = spatial derivatives of the evolved variables u.
 =    wrk = work space to for integrating the PDEs, used by a Runge-Kutta
 =          type integrator.
 = 
 =-----------------------------------------------------------------------=#
struct GridFields

    neqs :: Int64
    grid :: Grid
    u :: Array{Array{Float64, 1},1}
    u2 :: Array{Array{Float64, 1},1}
    dxu :: Array{Array{Float64, 1},1}
    wrk :: Array{Array{Float64, 1},1}
 
    function GridFields( neqs, grid )
   
        nx = grid.nx
        u = Array{Array{Float64, 1},1}(undef, neqs)
        u2 = Array{Array{Float64, 1},1}(undef, neqs)
        dxu = Array{Array{Float64, 1},1}(undef, neqs)
        wrk = Array{Array{Float64, 1},1}(undef, neqs)
        for i = 1:neqs
            u[i] = zeros(Float64,nx)
            u2[i] = zeros(Float64,nx)
            dxu[i] = zeros(Float64,nx)
            wrk[i] = zeros(Float64,nx)
        end
        new( neqs, grid, u, u2, dxu, wrk)

    end
end 

#=-----------------------------------------------------------------------
 =
 =  This routine computes the L2 norm of a function u.  
 =  This L2 Norm is an approximation to the integral over a uniform grid,
 =
 =                  [   Sum_i ( u_i^2 )  ]
 =  ||u|| =    Sqrt [------------------- ]
 =                  [       N - 1        ]
 =
 =-----------------------------------------------------------------------=#
function l2norm(u)
    s::Float64 = 0.0
    nx = length(u)

    for j = 1:nx
        s += u[j]*u[j]
    end
    return sqrt(s/(nx-1))
end

