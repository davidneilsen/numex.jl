using Printf
using WriteVTK
using ArgParse

include("Maxwell3D.jl")
import .Maxwell3D

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
    @printf("Step=%d, time=%g, |Ex|=%g, |Bz|=%g\n",0,time[1],Maxwell3D.l2norm(fields.u[1]), Maxwell3D.l2norm(fields.u[6]))
    vtk_grid(filename, fields.grid.x, fields.grid.y, fields.grid.z) do vtk
        vtk["Ex",VTKPointData()] = fields.u[1]
        vtk["Ey",VTKPointData()] = fields.u[2]
        vtk["Ez",VTKPointData()] = fields.u[3]
        vtk["Bx",VTKPointData()] = fields.u[4]
        vtk["By",VTKPointData()] = fields.u[5]
        vtk["Bz",VTKPointData()] = fields.u[6]
        vtk["time",VTKFieldData()] = time[1]
    end
    for i = 1:nt
        Maxwell3D.rk2_step!(Maxwell3D.maxwellEqs!, fields, time)
        if (mod(i,screenOutFreq)==0)
            @printf("Step=%d, time=%g, |Ex|=%g, |Bz|=%g\n",i,time[1],Maxwell3D.l2norm(fields.u[1]),Maxwell3D.l2norm(fields.u[6]))
        end
        if (mod(i,vtkOutFreq)==0)
            vtkFileCount += 1
            filename = @sprintf("maxwell_%05d",vtkFileCount)
            vtk_grid(filename, fields.grid.x, fields.grid.y, fields.grid.z) do vtk
                vtk["Ex",VTKPointData()] = fields.u[1]
                vtk["Ey",VTKPointData()] = fields.u[2]
                vtk["Ez",VTKPointData()] = fields.u[3]
                vtk["Bx",VTKPointData()] = fields.u[4]
                vtk["By",VTKPointData()] = fields.u[5]
                vtk["Bz",VTKPointData()] = fields.u[6]
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

    shp = [n, n, n]
    bbox = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]
    println("main:  nx = ",shp[1], " ny = ",shp[2], " nz = ",shp[3])
    grid = Maxwell3D.Grid(shp, bbox, cfl)
    println("grid:  nx = ",grid.shp[1], " ny = ",grid.shp[2], " nz = ",grid.shp[3])

    fields = Maxwell3D.GridFields(6, grid)
    Maxwell3D.init_data!(fields)
    
    evolve!(fields, nt, VTKOutFreq)
end

main()

