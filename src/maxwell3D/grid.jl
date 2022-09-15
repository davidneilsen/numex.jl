struct Grid

    bbox :: Array{Float64,1}
    shp  :: Array{Int64,1}

    x :: Vector{Float64}
    y :: Vector{Float64}
    z :: Vector{Float64}

    dx :: Float64
    dy :: Float64
    dz :: Float64
    dt :: Float64

    function Grid( sz, bbox, cfl)

        dx = (bbox[2]- bbox[1])/ (sz[1] - 1)
        dy = (bbox[4]- bbox[3])/ (sz[2] - 1)
        dz = (bbox[6]- bbox[5])/ (sz[3] - 1)
        c = 1.0
        dt = cfl / sqrt(1.0/dx^2+1.0/dy^2+1.0/dz^2) / c
        x = LinRange(bbox[1], bbox[2], sz[1])
        y = LinRange(bbox[3], bbox[4], sz[2])
        z = LinRange(bbox[5], bbox[6], sz[3])
        new( bbox, sz, x, y, z, dx, dy, dz, dt)

    end

end

struct GridFields

    neqs :: Int64
    grid :: Grid
    u :: Array{Array{Float64, 3},1}
    u2 :: Array{Array{Float64, 3},1}
    dxu :: Array{Array{Float64, 3},1}
    dyu :: Array{Array{Float64, 3},1}
    dzu :: Array{Array{Float64, 3},1}
    wrk :: Array{Array{Float64, 3},1}
 
    function GridFields( neqs, grid )
   
        nx, ny, nz = grid.shp[1], grid.shp[2], grid.shp[3]
        u = Array{Array{Float64, 3},1}(undef, neqs)
        u2 = Array{Array{Float64, 3},1}(undef, neqs)
        dxu = Array{Array{Float64, 3},1}(undef, neqs)
        dyu = Array{Array{Float64, 3},1}(undef, neqs)
        dzu = Array{Array{Float64, 3},1}(undef, neqs)
        wrk = Array{Array{Float64, 3},1}(undef, neqs)
        for i = 1:neqs
            u[i] = zeros(Float64,nx,ny,nz)
            u2[i] = zeros(Float64,nx,ny,nz)
            dxu[i] = zeros(Float64,nx,ny,nz)
            dyu[i] = zeros(Float64,nx,ny,nz)
            dzu[i] = zeros(Float64,nx,ny,nz)
            wrk[i] = zeros(Float64,nx,ny,nz)
        end
        new( neqs, grid, u, u2, dxu, dyu, dzu, wrk)

    end
end 

function l2norm(u::Array{Float64,3})
    s::Float64 = 0.0
    for j in CartesianIndices(u)
        s += u[j]*u[j]
    end
    nx, ny, nz = size(u)
    return sqrt(s/(nx*ny*nz))
end

