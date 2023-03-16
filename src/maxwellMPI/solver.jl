using Printf
using WriteVTK
using ArgParse
using MPI

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
    time = [0.0]
    vtkFileCount::Int64 = 0
    screenOutFreq = vtkOutFreq
    filename = @sprintf("maxwell_%05d",vtkFileCount)

    for i = 1:nt
        if (mod(i,screenOutFreq)==0)
            @printf("Step=%d, time=%g, |Bz|=%g\n",i,time[1],l2norm(fields.u[3]))
        end
    end

end

function main()
  
    root = 0

    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    numRanks = MPI.Comm_size(comm)

    print("@@ I am rank $(rank)) of $(numRanks))\n")

    n_params_int = 4
    n_params_double = 5
    params_int = Array{Int64}(undef, n_params_int)
    params_double = Array{Float64}(undef, n_params_double)

    if (rank == root) 
        println("@@ args = ",ARGS[1])
        fname = ARGS[1]
        println("@@ fname = ",fname)
        params_int, params_double = Maxwell2D.get_params(fname)
    end

    MPI.Barrier(comm)
    MPI.Bcast!(params_int, n_params_int, 0, comm)
    MPI.Bcast!(params_double, n_params_double, 0, comm)

    dims = MPI.Dims_create(numRanks, [0, 0])
    if (rank == root)
        println("@@ dims = ",dims)
    end

    nt = params_int[p_itmax]
    VTKOutFreq = params_int[p_vtkfreq]
    shp0 = [params_int[p_nx], params_int[p_ny]]

    bbox0 = Array{Float64,2}(undef, DIM, 2)
    bbox0[1,1] = params_double[p_xmin]
    bbox0[1,2] = params_double[p_xmax]
    bbox0[2,1] = params_double[p_ymin]
    bbox0[2,2] = params_double[p_ymax]
    cfl = params_double[p_cfl]
    ghostwidth = 3

    gh = Maxwell2D.GH(shp0, bbox0, cfl, numRanks, dims, ghostwidth, rank)

    # debug
    if (rank == 0)
        for i = 1:gh.boxcount
            @printf("min %d %d max %d %d dmin %g %g dmax %g %g, boxcount=%d\n",gh.ibox[i,1,1],gh.ibox[i,2,1],gh.ibox[i,1,2],gh.ibox[i,2,2],gh.cbox[i,1,1],gh.cbox[i,2,1],gh.cbox[i,1,2],gh.cbox[i,2,2],gh.boxcount)
        end
    end

    #grid = Maxwell2D.Grid(shp, bbox, cfl)
    #fields = Maxwell2D.GridFields(3, grid)
    #Maxwell2D.init_data!(fields)
    
    #evolve!(fields, nt, VTKOutFreq)
    MPI.Finalize()
end

main()

