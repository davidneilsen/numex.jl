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
    ut :: Array{Array{Float64, 1},1}
    dxu :: Array{Array{Float64, 1},1}
    k1 :: Array{Array{Float64, 1},1}
    k2 :: Array{Array{Float64, 1},1}
    k3 :: Array{Array{Float64, 1},1}
    k4 :: Array{Array{Float64, 1},1}

    A  :: Tridiagonal{Float64,Vector{Float64}}
    B  :: Array{Float64,2}
    F1 :: Tridiagonal{Float64,Vector{Float64}}
    F2 :: Array{Float64,2}
    order :: Int64
    method :: Int64
 
    function GridFields( neqs, grid, order, method=1 )
   
        nx = grid.nx
        u = Array{Array{Float64, 1},1}(undef, neqs)
        ut = Array{Array{Float64, 1},1}(undef, neqs)
        dxu = Array{Array{Float64, 1},1}(undef, neqs)
        k1 = Array{Array{Float64, 1},1}(undef, neqs)
        k2 = Array{Array{Float64, 1},1}(undef, neqs)
        k3 = Array{Array{Float64, 1},1}(undef, neqs)
        k4 = Array{Array{Float64, 1},1}(undef, neqs)
        for i = 1:neqs
            u[i] = zeros(Float64,nx)
            ut[i] = zeros(Float64,nx)
            dxu[i] = zeros(Float64,nx)
            k1[i] = zeros(Float64,nx)
            k2[i] = zeros(Float64,nx)
            k3[i] = zeros(Float64,nx)
            k4[i] = zeros(Float64,nx)
        end
        dx = grid.dx
        if order == 4
            A, B = dT4(nx, dx)
            F1, F2 = filterT6(nx, dx)
        elseif order == 6
            A, B = dT6(nx, dx)
            F1, F2 = filterT8(nx, dx)
        else
            println("Error.  Unknown order = ",order)
        end

        new( neqs, grid, u, ut, dxu, k1, k2, k3, k4, A, B, F1, F2, order, method)

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

