# The global grid
struct GH
    shp0   :: Array{Int64,1}    # global grid size
    boxcount :: Int64
    cbox0  :: Array{Float64,2}  # global coordinate bounding box

    ibox :: Array{Int64, 3}     # global indices for patches
    cbox :: Array{Float64, 3}   # global coordinate bounds for patches
    gibox :: Array{Int64, 3}    # global indices including ghostzones
    gcbox :: Array{Float64, 3}  # global coordinate bounds including ghostzones

    comm_partner :: Array{Int64, 1}
    length1 :: Array{Int64, 1}
    length2 :: Array{Int64, 1} 
    i1box :: Array{Int64, 3}
    i2box :: Array{Int64, 3}

    cfl  :: Float64
    dt   :: Float64
    dx0  :: Array{Float64, 1}   # dx in each dimension
    maxnum :: Int64             # number of mpi tasks
    ghostwidth :: Int64

    function GH( shp, bbox, cfl, size, dims, ghostwidth,rank)
        dx0 = Array{Float64,1}(undef, DIM)
        for i = 1:DIM
            dx0[i] = (bbox[i,2] - bbox[i,1])/(shp[i]-1)
        end

        dt = dx0[1] * cfl

        ibox = zeros(Int64,size,DIM,2)
        gibox = zeros(Int64,size,DIM,2)

        boxcount = 0
    
        minx = 1
        maxx = shp[1]
        miny = 1
        maxy = shp[2]
    
        lengthi = round(Int,(maxx - minx)/dims[1])
        lengthj = round(Int,(maxy - miny)/dims[2])
    
        numx = dims[1]
        numy = dims[2]
    
        if ( (numx > 1 || numy > 1) && (size >= (boxcount + numx*numy)) )
            for j = 1:numy
                for i = 1:numx
                    boxcount += 1
    
                    ibox[boxcount,1,1] = minx + (i-1)*lengthi
                    ibox[boxcount,1,2] = minx + (i)*lengthi
                    ibox[boxcount,2,1] = miny + (j-1)*lengthj
                    ibox[boxcount,2,2] = miny + (j)*lengthj
    
                    gibox[boxcount,1,1] = minx + (i-1)*lengthi
                    gibox[boxcount,1,2] = minx + (i)*lengthi
                    gibox[boxcount,2,1] = miny + (j-1)*lengthj
                    gibox[boxcount,2,2] = miny + (j)*lengthj
    
                    if i == 1 ibox[boxcount,1,1] = minx end
                    if j == 1 ibox[boxcount,2,1] = miny end
                    if i == numx ibox[boxcount,1,2] = maxx end
                    if j == numy ibox[boxcount,2,2] = maxy end
    
                    if i == 1 gibox[boxcount,1,1] = minx end
                    if j == 1 gibox[boxcount,2,1] = miny end
                    if i == numx gibox[boxcount,1,2] = maxx end
                    if j == numy gibox[boxcount,2,2] = maxy end
    
                    # Add ghostzones
                    if ibox[boxcount,1,1] != minx ibox[boxcount,1,1] -= ghostwidth end
                    if ibox[boxcount,1,2] != maxx ibox[boxcount,1,2] += ghostwidth end
                    if ibox[boxcount,2,1] != miny ibox[boxcount,2,1] -= ghostwidth end
                    if ibox[boxcount,2,2] != maxy ibox[boxcount,2,2] += ghostwidth end

                    if rank==0
                    @printf("i,j=(%d,%d), numx/y=%d,%d, size=%d\n",i,j,numx,numy,size)
                    @printf("lengthi=%d, lengthj=%d\n",lengthi,lengthj)
                    println(ibox)
                    end
    
                end
            end
        else
            boxcount += 1
    
            ibox[boxcount,1,1] = minx
            ibox[boxcount,1,2] = maxx
            ibox[boxcount,2,1] = miny
            ibox[boxcount,2,2] = maxy
            
            gibox[boxcount,1,1] = minx
            gibox[boxcount,1,2] = maxx
            gibox[boxcount,2,1] = miny
            gibox[boxcount,2,2] = maxy
        end
    
        cbox = zeros(Float64,boxcount,DIM,2)
        gcbox = zeros(Float64,boxcount,DIM,2)
        comm_partner = zeros(Int64,boxcount)
        length1 = zeros(Int64,boxcount)
        length2 = zeros(Int64,boxcount)
        i1box = zeros(Int64,boxcount,DIM,2)
        i2box = zeros(Int64,boxcount,DIM,2)

        for i = 1:boxcount
            for j = 1:DIM
                cbox[i,j,1] = bbox[j,1] + dx0[j]*(ibox[i,j,1]-1)
                cbox[i,j,2] = bbox[j,1] + dx0[j]*(ibox[i,j,2]-1)
                gcbox[i,j,1] = bbox[j,1] + dx0[j]*(gibox[i,j,1]-1)
                gcbox[i,j,2] = bbox[j,1] + dx0[j]*(gibox[i,j,2]-1)
            end
        end
         
        new( shp, boxcount, bbox, ibox, cbox, gibox, gcbox, comm_partner, length1, length2, i1box, i2box, cfl, dt, dx0, size, ghostwidth)
    end
end

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

    function Grid( shp, bbox, cfl)

        nx = shp[1]
        ny = shp[2]
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

function l2norm_global(u::Array{Float64,2})
    s::Float64 = 0.0
    nx = gh.nx
    ny = gh.ny
    ng = gh.ghostwidth
    hx = gh.hx
    hy = gh.hy

    for j = ng:ny-ng
        for i = ng:nx-ng
            s += u[j]*u[j]
        end
    end
    s *= hx*hy 

end