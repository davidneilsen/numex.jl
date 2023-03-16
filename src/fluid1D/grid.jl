struct Pars

end


struct Grid

    bbox :: Vector{Float64}
    nx   :: Int64

    x   :: Vector{Float64}
    dx  :: Float64
    dt  :: Float64

    function create_coords(xmin, xmax, sz)
       x = zeros(Float64, sz)
       dx = (xmax - xmin)/(sz -1)
       for i = 1:sz
           x[i] = xmin + (i-1)*dx
       end
       return x, dx
    end

    function Grid( nx, bbox, cfl)
        x, dx = create_coords(bbox[1], bbox[2], nx)

        # set dt
        c::Float64 = 1.0
        tmp = 1.0/dx^2
        dt = cfl / sqrt(sum(tmp)) / c
        #println("dt = ",dt," cfl = ",cfl," tmp = ",tmp)

        new( bbox, nx, x, dx, dt)

    end

end
 
struct GridFields

    neqs :: Int64
    grid :: Grid
    u :: Array{Vector{Float64},1}
    v :: Array{Vector{Float64},1}
    u2 :: Array{Vector{Float64},1}
    du :: Array{Vector{Float64},1}
    wrk :: Array{Array{Vector{Float64},1},}
 
    function GridFields( neqs, grid, nstage )

        nx = grid.nx
        u = Array{Vector{Float64}}(undef, neqs)
        v = Array{Vector{Float64}}(undef, neqs)
        u2 = Array{Vector{Float64}}(undef, neqs)
        du = Array{Vector{Float64}}(undef, neqs)
        wrk = Array{Array{Vector{Float64}},1}(undef, nstage)
        pars = Dict( "gamma" => 1.5, "idst_rhoL" => 1.0, "idst_vxL" => 0.0, "idst_PL" => 1000.0, 
        "idst_rhoR" => 1.0, "idst_vxR" => 0.0, "idst_PR" => 0.01 )

        println("nx = ",nx)
        for i = 1:nstage
            wrk[i] = Array{Vector{Float64},1}(undef, neqs)
        end

        for i = 1:neqs
            u[i] = Vector{Float64}(undef,nx)
            v[i] = Vector{Float64}(undef,nx)
            u2[i] = Vector{Float64}(undef,nx)
            du[i] = Vector{Float64}(undef,nx)
            for j = 1:nstage
               wrk[j][i] = Vector{Float64}(undef,nx)
            end
        end
        new( neqs, grid, u, v, u2, du, wrk)

    end
end 

function l2norm(u::Vector{Float64})
    s::Float64 = 0.0
    nx = length(u)
    for j in 1:nx
        s += u[j]*u[j]
    end
    return sqrt(s/j)
end

