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

function l2norm(u)
    s::Float64 = 0.0
    nx = length(u)

    for j = 1:nx
        s += u[j]*u[j]
    end
    return sqrt(s/nx)
end

