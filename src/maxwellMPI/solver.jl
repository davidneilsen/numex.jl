using Printf
using WriteVTK
using ArgParse
using MPI
using BlockBandedMatrices, BandedMatrices, ArrayLayouts, FillArrays, LazyBandedMatrices, LazyArrays

include("pars.jl")
include("Maxwell2D.jl")
using .Maxwell2D

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx", "-n"
            help = "number of points in each direction"
            arg_type = Int
            default = 201
       "--vtkfreq", "-f"
            help = "Frequency of VTK output"
            arg_type = Int
            default = 5
         "--cfl", "-c"
            help = "Courant number"
            arg_type = Float64
            default = 0.1
        "nsteps"
            help = "Number of steps"
            arg_type = Int
            required = true
 
    end
    return parse_args(s)
end

function evolve!(fields, nt, vtkOutFreq)
    time = zeros(Float64,1)
    vtkFileCount::Int64 = 0
    screenOutFreq = vtkOutFreq

    gh = fields.gh
    lrank = gh.gridID
    numRanks = gh.numRanks
    lx = gh.lcoords[1]
    ly = gh.lcoords[2]

    extent = Vector{Tuple{UnitRange{Int64}, UnitRange{Int64}}}(undef,gh.boxcount)
    for k = 1:gh.boxcount
        extent[k] =(gh.gibox[k,1,1]:gh.gibox[k,1,2], gh.gibox[k,2,1]:gh.gibox[k,2,2]) 
    end
    #print("extents is ", extent)

    gnrm = Maxwell2D.l2norm_global(fields.u[1],gh)
    if lrank == 1
        @printf("Step=%d, time=%g, |Ex|=%g\n",0, time[1], gnrm)
    end
    filename = @sprintf("maxwell_%05d",vtkFileCount)
    pvtk_grid(filename, lx, ly; part=lrank, nparts=numRanks, extents=extent) do pvtk
        pvtk["proc",VTKPointData()] = fields.proc
        pvtk["Ex",VTKPointData()] = fields.u[1]
        pvtk["Ey",VTKPointData()] = fields.u[2]
        pvtk["Hz",VTKPointData()] = fields.u[3]
        pvtk["time",VTKFieldData()] = time[1]
    end

    for i = 1:nt
        #@printf("@@@ Step=%d, time=%g, |Bz|=%g\n",i,time[1],l2norm(fields.u[3]))
        Maxwell2D.rk2_step!(Maxwell2D.maxwell_TE!, fields, time)
        gnrm = Maxwell2D.l2norm_global(fields.u[1],gh)
        if lrank == 1 && mod(i,screenOutFreq)==0
            @printf("Step=%d, time=%g, |Ex|=%g\n",i,time[1],gnrm)
        end
        if (mod(i,vtkOutFreq)==0)
            vtkFileCount += 1
            filename = @sprintf("maxwell_%05d",vtkFileCount)
            pvtk_grid(filename, lx, ly; part=lrank, nparts=numRanks, extents=extent) do pvtk
                pvtk["proc",VTKPointData()] = fields.proc
                pvtk["Ex",VTKPointData()] = fields.u[1]
                pvtk["Ey",VTKPointData()] = fields.u[2]
                pvtk["Hz",VTKPointData()] = fields.u[3]
                pvtk["time",VTKFieldData()] = time[1]
            end
        end
    end

end

function main()
    ltrace = false
  
    root = 0

    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    numRanks = MPI.Comm_size(comm)

    # print("@@ I am rank $(rank)) of $(numRanks))\n")

    n_params_int = 6
    n_params_double = 5
    params_int = Array{Int64}(undef, n_params_int)
    params_double = Array{Float64}(undef, n_params_double)

    if (rank == root) 
        # println("@@ args = ",ARGS[1])
        fname = ARGS[1]
        # println("@@ fname = ",fname)
        params_int, params_double = Maxwell2D.get_params(fname)
    end

    MPI.Barrier(comm)
    MPI.Bcast!(params_int, n_params_int, 0, comm)
    MPI.Bcast!(params_double, n_params_double, 0, comm)

    dims = MPI.Dims_create(numRanks, [0, 0])
    if (ltrace && rank == root)
        println("@@ dims = ",dims)
    end

    nt = params_int[p_itmax]
    VTKOutFreq = params_int[p_vtkfreq]
    shp0 = [params_int[p_nx], params_int[p_ny]]
    dtype = params_int[p_dtype]

    bbox0 = Array{Float64,2}(undef, DIM, 2)
    bbox0[1,1] = params_double[p_xmin]
    bbox0[1,2] = params_double[p_xmax]
    bbox0[2,1] = params_double[p_ymin]
    bbox0[2,2] = params_double[p_ymax]
    cfl = params_double[p_cfl]
    ghostwidth = params_int[p_ghostwidth]

    gh = Maxwell2D.GH(shp0, bbox0, cfl, numRanks, dims, ghostwidth, comm, rank, dtype)

    if numRanks > 1 
        Maxwell2D.set_communication(gh) 
    end

    # debug
    if (ltrace && rank == 0)
        for i = 1:gh.boxcount
            #@printf("min %d %d max %d %d dmin %g %g dmax %g %g, boxcount=%d\n",gh.ibox[i,1,1],gh.ibox[i,2,1],gh.ibox[i,1,2],gh.ibox[i,2,2],gh.cbox[i,1,1],gh.cbox[i,2,1],gh.cbox[i,1,2],gh.cbox[i,2,2],gh.boxcount)
            @printf("gridID %d: ibox (%d, %d), (%d, %d) | ngibox (%d, %d), cbox(%g, %g)\n",i, gh.ibox[i,1,1], gh.ibox[i,1,2], gh.ibox[i,2,1], gh.ibox[i,2,2], gh.gibox[i,1,1],gh.gibox[i,1,2],gh.cbox[i,1,1],gh.cbox[i,1,2])
        end
    end

    # MPI.Barrier(comm)

    if ltrace
        for k = 1:gh.boxcount
            if rank==(k-1)
                @printf("gridID %d: lshp = %d, %d, | bbox = (%d, %d), (%d, %d) | x = (%g, %g), y = (%g, %g)\n",k, gh.lshp[1],gh.lshp[2],gh.lcbox[1,1],gh.lcbox[1,2],gh.lcbox[2,1],gh.lcbox[2,2],gh.lcoords[1][1],gh.lcoords[1][gh.lshp[1]],gh.lcoords[2][1],gh.lcoords[2][gh.lshp[2]])
                MPI.Barrier(comm)
                for i = 1:gh.comm_count[1]
                    @printf("gridID %d: from %d recv (%d, %d), (%d, %d) len = %d | send (%d, %d), (%d, %d)len = %d\n",k, gh.comm_partner[i], gh.irecbox[i,1,1], gh.irecbox[i,1,2], gh.irecbox[i,2,1], gh.irecbox[i,2,2], gh.lensend[i], gh.isendbox[i,1,1],gh.isendbox[i,1,2],gh.isendbox[i,2,1],gh.isendbox[i,2,2],gh.lensend[i])
                end
            end
        end
    end


    fields = Maxwell2D.GridFields(3, gh)
    Maxwell2D.waveguide_init_data!(fields)
    Maxwell2D.grid_sync!(fields.u, gh, comm)
    
    evolve!(fields, nt, VTKOutFreq)
    MPI.Finalize()
end

main()

