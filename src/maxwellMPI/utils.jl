using MPI

function nint(x)
    if (x >= 0.0)

        n = round(Int, x)
    else
        n = round(Int, x)
    end
    return n
end

function set_communication(gh::GH)

    gh.comm_count[1] = 0
    for j = 1:gh.boxcount
        if j != gh.gridID
            if grid_intersect(gh.gridID,j,gh.gibox) == 1
                gh.comm_count[1] += 1
                gh.comm_partner[gh.comm_count[1]] = j
                #@printf("Rank %d: partner %d\n",gh.gridID,j)
            end
        end
    end

    length = ones(Int64,1)
    icombox = Array{Int64,2}(undef,2,2)
    for i = 1:gh.comm_count[1]
        find_intersection!(icombox, length, gh.gridID, gh.comm_partner[i], gh.ibox, gh.gibox)
        for j = 1:2
            for k = 1:2
                gh.irecbox[i,j,k] = icombox[j,k]
            end
        end
        gh.lenrec[i] = length[1]
        find_intersection!(icombox, length, gh.comm_partner[i], gh.gridID, gh.ibox, gh.gibox)
        for j = 1:2
            for k = 1:2
                gh.isendbox[i,j,k] = icombox[j,k]
            end
        end
        gh.lensend[i] = length[1]
    end

end

function grid_intersect(g1, g2, gibox)

    LTRACE = false
    LTRACE2 = false 

    if LTRACE
        @printf("grid_intersect: checking grids %d and %d\n",g1,g2)
    end

    if LTRACE2
        @printf("... g1=%d,k=1,gibox(g1,min/max)=%d/%d, g2=%d,k=1,gibox(g2,min/max)=%d/%d",g1,gibox[g1,1,1],gibox[g1,1,2],g2,gibox[g2,1,1],gicbox[g2,1,2])
        @printf("... g1=%d,k=2,gibox(g1,min/max)=%d/%d, g2=%d,k=2,gibox(g2,min/max)=%d/%d",g1,gibox[g1,1,1],gibox[g1,2,2],g2,gicbox[g2,1,1],gicbox[g2,2,2])
    end

    for k = 1:DIM
        if (gibox[g1,k,1] >= gibox[g2,k,1])
            if (gibox[g1,k,1] > gibox[g2,k,2])
                # no intersection
                if LTRACE2 
                    @printf("...grid_intersect: grids(%d, %d) no intersection 1, k=%d\n",g1,g2,k) 
                    #@printf("......k=%d, g1 min: %d, g2 max: %d\n",k,gibox[g1,k,1],gibox[g2,k,2])
                end
                return 0
            end
        else
            if (gibox[g1,k,2] < gibox[g2,k,1])
                # no intersection
                if LTRACE2
                    @printf("...grid_intersect: grids(%d, %d) no intersection 2, k=%d\n",g1,g2,k) 
                    #@printf("......k=%d, g1 max: %d, g2 min: %d\n",k,gibox[g1,k,2],gibox[g2,k,1])
                end
                return 0
            end
        end
    end

    if LTRACE 
        @printf("...grid_intersect: grids(%d, %d) intersection!\n",g1,g2) 
    end

    return 1
end

function find_intersection!(icombox, length, g1, g2, ibox, gibox)
    LTRACE = false
    LTRACE2 = false

    if LTRACE @printf("find_intersection!\n") end

#    if LTRACE2
#        @printf("g1 %d, g2 %d, rank %d, ngdminx %g, ngdmaxx %g\n",g1,g2);
#    end

    for k = 1:DIM
        if gibox[g1,k,1] >= ibox[g2,k,1]
            if gibox[g1,k,1] > ibox[g2,k,2] return 0 end
            icombox[k,1] = gibox[g1,k,1]
            icombox[k,2] = ibox[g2,k,2]
        else
            if gibox[g1,k,2] < ibox[g2,k,1] return 0 end
            icombox[k,1] = ibox[g2,k,1]
            icombox[k,2] = gibox[g1,k,2]
        end
    end

    length[1] = (icombox[1,2] - icombox[1,1] + 1)*(icombox[2,2] - icombox[2,1] + 1)
    if (LTRACE)
        @printf("g1 %d, g2 %d, icombox: %d, %d, %d, %d\n",g1,g2,icombox[1,1],icombox[1,2],icombox[2,1],icombox[2,2])
    end
end

function grid_sync!(u, gh, comm)

    neqs = length(u)
    if gh.numRanks == 1 return 0 end

    comm_count = gh.comm_count[1]

    send_buffer = Array{Array{Float64,1},1}(undef, comm_count)
    recv_buffer = Array{Array{Float64,1},1}(undef, comm_count)

    for i = 1:comm_count
        send_buffer[i] = Array{Float64,1}(undef, neqs*gh.lensend[i])
        recv_buffer[i] = Array{Float64,1}(undef, neqs*gh.lenrec[i])
    end

    sreq = Vector{MPI.Request}(undef, comm_count)
    rreq = Vector{MPI.Request}(undef, comm_count)

    for k = 1:comm_count
        cpk = gh.comm_partner[k] - 1

        #@printf("rank: %d >>> send to %d. len=%d, (%d, %d) (%d, %d)\n",gh.rank,cpk,gh.lensend[k],gh.isendbox[k,1,1],gh.isendbox[k,1,2],gh.isendbox[k,2,1],gh.isendbox[k,2,2])

        rreq[k] = MPI.Irecv!(recv_buffer[k], cpk, cpk, comm)

        ls = length(send_buffer[k])
        for m = 1:neqs
            n = 0
            for j = gh.isendbox[k,2,1]:gh.isendbox[k,2,2]
                jj = j - gh.loffset[2]
                for i = gh.isendbox[k,1,1]:gh.isendbox[k,1,2]
                    ii = i - gh.loffset[1]
                    n += 1
                    if (lii=n+(m-1)*gh.lensend[k])>ls
                        @printf(">>> problem! ls=%d,lii=%d, n=%d, k=%d,m=%d,i=%d,j=%d,lsend=%d,(%d,%d),(%d,%d)\n",ls,lii,n,k,m,i,j,gh.lensend[k],gh.isendbox[k,1,1],gh.isendbox[k,1,2],gh.isendbox[k,2,1],gh.isendbox[k,2,2])
                    end
                    send_buffer[k][n+(m-1)*gh.lensend[k]] = u[m][ii,jj]
                end
            end
        end
        sreq[k] = MPI.Isend(send_buffer[k], cpk, gh.rank, comm)
        #@printf("rank: %d >>> send to %d. len=%d, (%d, %d) (%d, %d), tag=%d, |u|=%g\n",gh.rank,cpk,gh.lensend[k],gh.isendbox[k,1,1],gh.isendbox[k,1,2],gh.isendbox[k,2,1],gh.isendbox[k,2,2],gh.rank,l2norm(send_buffer[k]))
    end
                    
    MPI.Waitall(sreq)
    MPI.Waitall(rreq)

    for k = 1:comm_count
        cpk = gh.comm_partner[k] - 1
        #@printf("rank: %d >>> received from %d. len=%d, (%d, %d) (%d, %d), tag=%d, |u|=%g\n",gh.rank,cpk,gh.lenrec[k],gh.irecbox[k,1,1],gh.irecbox[k,1,2],gh.irecbox[k,2,1],gh.irecbox[k,2,2],gh.rank,l2norm(recv_buffer[k]))
        for m = 1:neqs
            n = 0
            for j = gh.irecbox[k,2,1]:gh.irecbox[k,2,2]
                jj = j - gh.loffset[2]
                for i = gh.irecbox[k,1,1]:gh.irecbox[k,1,2]
                    ii = i - gh.loffset[1]
                    n += 1
                    u[m][ii,jj] = recv_buffer[k][n+(m-1)gh.lenrec[k]]
                end
            end
        end
    end
    
end

function l2norm(u::Vector{Float64})
    s::Float64 = 0.0
    nx = length(u)
    for j = 1:nx
        s += u[j]*u[j]
    end
    return sqrt(s/nx)
end

function l2norm(u::Array{Float64,2})
    s::Float64 = 0.0
    nx, ny = size(u)
    for j = 1:ny
        for i = 1:nx
            s += u[j]*u[j]
        end
    end
    return sqrt(s/(nx*ny))
end

function l2norm_global(u::Array{Float64,2}, gh)

    LTRACE = false
    s::Float64 = 0.0

    gridID = gh.gridID
    jb = gh.ibox[gridID,2,1] - gh.loffset[2]
    je = gh.ibox[gridID,2,2] - gh.loffset[2]
    ib = gh.ibox[gridID,1,1] - gh.loffset[1]
    ie = gh.ibox[gridID,1,2] - gh.loffset[1]

    hx = gh.dx0[1]
    hy = gh.dx0[2]

    if LTRACE
        @printf("proc=%d, irange=(%d, %d), jrange=(%d, %d), ibox=(%d,%d),(%d,%d), gibox=(%d,%d),(%d,%d)\n",gh.rank,ib,ie,jb,je,gh.ibox[gridID,1,1],gh.ibox[gridID,1,2],gh.ibox[gridID,2,1],gh.ibox[gridID,2,2],gh.gibox[gridID,1,1],gh.gibox[gridID,1,2],gh.gibox[gridID,2,1],gh.gibox[gridID,2,2])
    end
 
    for j = jb:je
        for i = ib:ie
            s += u[i,j]*u[i,j]
        end
    end
    s *= hx*hy 

    gs = MPI.Reduce(s, MPI.SUM, gh.comm; root=0)

    MPI.Barrier(gh.comm)
    if gh.rank == 0
        return sqrt(gs)
    else
        return -1.0
    end

end
