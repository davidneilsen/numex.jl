using Printf
using WriteVTK
using ArgParse

include("Maxwell.jl")
import .Maxwell

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx", "-n"
            help = "number of points in each direction"
            arg_type = Int
            default = 401
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
    @printf("Step=%d, time=%g, |Bz|=%g\n",0,time[1],Maxwell.l2norm(fields.u[3]))
    vtk_grid(filename, fields.grid.x, fields.grid.y) do vtk
        vtk["Ex",VTKPointData()] = fields.u[1]
        vtk["Ey",VTKPointData()] = fields.u[2]
        vtk["Hz",VTKPointData()] = fields.u[3]
        vtk["time",VTKFieldData()] = time[1]
    end
    for i = 1:nt
        Maxwell.rk2_step!(Maxwell.maxwell_TE!, fields, time)
        if (mod(i,screenOutFreq)==0)
            @printf("Step=%d, time=%g, |Bz|=%g\n",i,time[1],Maxwell.l2norm(fields.u[3]))
        end
        if (mod(i,vtkOutFreq)==0)
            vtkFileCount += 1
            filename = @sprintf("maxwell_%05d",vtkFileCount)
            vtk_grid(filename, fields.grid.x, fields.grid.y) do vtk
                vtk["Ex",VTKPointData()] = fields.u[1]
                vtk["Ey",VTKPointData()] = fields.u[2]
                vtk["Hz",VTKPointData()] = fields.u[3]
                vtk["time",VTKFieldData()] = time[1]
            end
        end
    end
end

function main()
  
    params = parse_commandline()
    n = params["nx"]
    nt = params["nsteps"]
    cfl = params["cfl"]
    VTKOutFreq = params["vtkfreq"]

    nx = ny = n
    bbox = [-10.0, 10.0, -10.0, 10.0]
    println("main:  nx = ",nx, " ny = ",ny)
    grid = Maxwell.Grid(nx, ny, bbox, cfl)
    println("main grid:  nx = ",grid.nx, " ny = ",grid.ny)

    fields = Maxwell.GridFields(3, grid)
    Maxwell.init_data!(fields)
    
    evolve!(fields, nt, VTKOutFreq)
end

main()

