using MPI
using BlockBandedMatrices, BandedMatrices, ArrayLayouts, FillArrays, LazyBandedMatrices, LazyArrays

# The global grid
struct GH
    shp0   :: Array{Int64,1}    # global grid size
    boxcount :: Int64
    cbox0  :: Array{Float64,2}  # global coordinate bounding box

    ibox :: Array{Int64, 3}     # global indices for patches
    cbox :: Array{Float64, 3}   # global coordinate bounds for patches
    gibox :: Array{Int64, 3}    # global indices including ghostzones
    gcbox :: Array{Float64, 3}  # global coordinate bounds including ghostzones

    comm :: MPI.Comm
    comm_count :: Vector{Int64}
    comm_partner :: Array{Int64, 1}
    lensend :: Array{Int64, 1}
    lenrec :: Array{Int64, 1} 
    isendbox :: Array{Int64, 3}
    irecbox :: Array{Int64, 3}

    cfl  :: Float64
    dt   :: Float64
    dx0  :: Array{Float64, 1}   # dx in each dimension
    numRanks :: Int64             # number of mpi tasks
    ghostwidth :: Int64
    rank :: Int64
    gridID :: Int64

    lshp :: Vector{Int64}
    loffset :: Vector{Int64}
    lcoords :: Vector{Vector{Float64}}
    lcbox :: Array{Float64,DIM}

    dtype :: Int64
    dissipation :: Int64

    function GH( shp, bbox, cfl, size, dims, ghostwidth, comm, rank, difftype, dissipation)

        LTRACE = false

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
                    ibox[boxcount,1,2] = minx + (i)*lengthi - 1
                    ibox[boxcount,2,1] = miny + (j-1)*lengthj
                    ibox[boxcount,2,2] = miny + (j)*lengthj - 1
    
                    gibox[boxcount,1,1] = minx + (i-1)*lengthi
                    gibox[boxcount,1,2] = minx + (i)*lengthi - 1
                    gibox[boxcount,2,1] = miny + (j-1)*lengthj
                    gibox[boxcount,2,2] = miny + (j)*lengthj - 1
    
                    if i == 1 ibox[boxcount,1,1] = minx end
                    if j == 1 ibox[boxcount,2,1] = miny end
                    if i == numx ibox[boxcount,1,2] = maxx end
                    if j == numy ibox[boxcount,2,2] = maxy end
    
                    if i == 1 gibox[boxcount,1,1] = minx end
                    if j == 1 gibox[boxcount,2,1] = miny end
                    if i == numx gibox[boxcount,1,2] = maxx end
                    if j == numy gibox[boxcount,2,2] = maxy end
    
                    # Add ghostzones
                    if ibox[boxcount,1,1] != minx gibox[boxcount,1,1] -= ghostwidth end
                    if ibox[boxcount,1,2] != maxx gibox[boxcount,1,2] += ghostwidth end
                    if ibox[boxcount,2,1] != miny gibox[boxcount,2,1] -= ghostwidth end
                    if ibox[boxcount,2,2] != maxy gibox[boxcount,2,2] += ghostwidth end

                    if LTRACE && rank==0
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
        comm_count = zeros(Int64,2)
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

        # gridID of the local proc
        gridID = rank+1
         
        # Set local grid size
        lshp = zeros(Int64,DIM)
        loffset = zeros(Int64,DIM)
        lcoords = Vector{Vector{Float64}}(undef,DIM)
        lcbox = zeros(Float64,DIM,2)

        for k = 1:DIM
            lshp[k] = gibox[gridID,k,2] - gibox[gridID,k,1] + 1
            loffset[k] = gibox[gridID,k,1] - 1
            lcbox[k,1] = gcbox[gridID,k,1]
            lcbox[k,2] = gcbox[gridID,k,2]
            lcoords[k] = LinRange(lcbox[k,1], lcbox[k,2], lshp[k])
        end

        dtype = difftype
        diss = dissipation

        new( shp, boxcount, bbox, ibox, cbox, gibox, gcbox, comm, comm_count, comm_partner, length1, length2, i1box, i2box, cfl, dt, dx0, size, ghostwidth, rank, gridID, lshp, loffset, lcoords, lcbox, dtype, diss)
    end
end

struct DerivVars

    dtype :: Int
    Dx1d :: Matrix
    Dy1d :: Matrix
    Dx2d :: BandedBlockBandedMatrix
    Dy2d :: BandedBlockBandedMatrix
    u1x :: Vector{Float64}
    u1y :: Vector{Float64}
    du1x :: Vector{Float64}
    du1y :: Vector{Float64}

    Fx1d :: BandedMatrix
    Fy1d :: BandedMatrix

    dissipation :: Int

    function DerivVars(shp, dx0, dtype, diss)
        nx = shp[1]
        ny = shp[2]
        hx = dx0[1]
        hy = dx0[2]

        Px, Qx = derivKim4(nx,hx)
        Py, Qy = derivKim4(ny,hy)
        Dx = Px \ Qx
        Dy = Py \ Qy
        Dx2d = BandedBlockBandedMatrix(Kron(Eye(ny), Dx))
        Dy2d = BandedBlockBandedMatrix(Kron(Dy, Eye(nx)))

        Ax, Bx = filterKim6(nx)
        Ay, By = filterKim6(ny)
        Fx = Ax \ Bx
        Fy = Ay \ By

        u1x = Vector{Float64}(undef, nx)
        u1y = Vector{Float64}(undef, ny)
        du1x = Vector{Float64}(undef, nx)
        du1y = Vector{Float64}(undef, ny)

        new(dtype, Dx, Dy, Dx2d, Dy2d, u1x, u1y, du1x, du1y, Fx, Fy, diss)
    end 
end

struct GridFields

    neqs :: Int64
    gh :: GH
    u :: Array{Array{Float64, 2},1}
    u2 :: Array{Array{Float64, 2},1}
    dxu :: Array{Array{Float64, 2},1}
    dyu :: Array{Array{Float64, 2},1}
    wrk :: Array{Array{Array{Float64, 2},1},1}
    proc :: Array{Float64,2}
    dvars :: DerivVars
 
    function GridFields( neqs, gh::GH)
   
        nx = gh.lshp[1]
        ny = gh.lshp[2]
        u = Array{Array{Float64, 2},1}(undef, neqs)
        u2 = Array{Array{Float64, 2},1}(undef, neqs)
        dxu = Array{Array{Float64, 2},1}(undef, neqs)
        dyu = Array{Array{Float64, 2},1}(undef, neqs)

        for i = 1:neqs
            u[i] = zeros(Float64,nx,ny)
            u2[i] = zeros(Float64,nx,ny)
            dxu[i] = zeros(Float64,nx,ny)
            dyu[i] = zeros(Float64,nx,ny)
        end

        nrk = 4
        wrk = Array{Array{Array{Float64, 2},1},1}(undef,nrk)
        for m = 1:nrk
            wrk[m] = Array{Array{Float64, 2},1}(undef, neqs)
        end
        for m = 1:nrk, i = 1:neqs
            wrk[m][i] = zeros(Float64,nx,ny)
        end

        proc = zeros(nx,ny)
        dvars = DerivVars(gh.lshp, gh.dx0, gh.dtype, gh.dissipation)

        new( neqs, gh, u, u2, dxu, dyu, wrk, proc, dvars)

    end
end 
