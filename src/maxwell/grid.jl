struct Grid

    xmin :: Float64
    xmax :: Float64
    ymin :: Float64
    ymax :: Float64

    nx :: Int64
    ny :: Int64

    x :: Vector{Float64}
    y :: Vector{Float64}

    dx :: Float64
    dy :: Float64
    dt :: Float64

    function Grid( nx0, ny0, bbox, cfl)

        nx = nx0
        ny = ny0
        xmin = bbox[1]
        xmax = bbox[2]
        ymin = bbox[3]
        ymax = bbox[4]
        dx = (xmax - xmin)/ (nx - 1)
        dy = (ymax - ymin)/ (ny - 1)
        c = 1.0
        dt = cfl / sqrt(1/dx^2+1/dy^2) / c
        x = LinRange(xmin, xmax, nx)
        y = LinRange(ymin, ymax, ny)
        new( xmin, xmax, ymin, ymax, nx, ny, x, y, dx, dy, dt)

    end

end

struct GridFields

    neqs :: Int64
    grid :: Grid
    u :: Array{Array{Float64, 2},1}
    u2 :: Array{Array{Float64, 2},1}
    dxu :: Array{Array{Float64, 2},1}
    dyu :: Array{Array{Float64, 2},1}
    wrk :: Array{Array{Float64, 2},1}
 
    function GridFields( neqs, grid )
   
        nx, ny = grid.nx, grid.ny
        u = Array{Array{Float64, 2},1}(undef, neqs)
        u2 = Array{Array{Float64, 2},1}(undef, neqs)
        dxu = Array{Array{Float64, 2},1}(undef, neqs)
        dyu = Array{Array{Float64, 2},1}(undef, neqs)
        wrk = Array{Array{Float64, 2},1}(undef, neqs)
        for i = 1:neqs
            u[i] = zeros(Float64,nx,ny)
            u2[i] = zeros(Float64,nx,ny)
            dxu[i] = zeros(Float64,nx,ny)
            dyu[i] = zeros(Float64,nx,ny)
            wrk[i] = zeros(Float64,nx,ny)
        end
        new( neqs, grid, u, u2, dxu, dyu, wrk)

    end
end 

function l2norm(u::Array{Float64,2})
    s::Float64 = 0.0
    for j in CartesianIndices(u)
        s += u[j]*u[j]
    end
    nx, ny =size(u)
    return sqrt(s/(nx*ny))
end

