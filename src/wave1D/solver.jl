using Printf
using ArgParse
using Plots

include("Wave1D.jl")
import .Wave1D

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx", "-n"
            help = "number of points in each direction"
            arg_type = Int
            default = 101
       "--outfreq", "-f"
            help = "Frequency of output"
            arg_type = Int
            default = 2
         "--cfl", "-c"
            help = "Courant number"
            arg_type = Float64
            default = 0.5
        "nsteps"
            help = "Number of steps"
            arg_type = Int
            required = true
 
    end
    return parse_args(s)
end

function evolve!(fields, nt, outFreq)
    time = [0.0]
    api = Animation()
    aphi = Animation()
    xrange = (fields.grid.x[1],fields.grid.x[end])
    yrange = (-0.5,0.5)

    screenOutFreq = outFreq
    @printf("Step=%d, time=%g, |Pi|=%g\n",0.0,time[1],Wave1D.l2norm(fields.u[1]))
    for i = 1:nt
        Wave1D.rk2_step!(Wave1D.waveEqs!, fields, time)
        if (mod(i,screenOutFreq)==0)
            @printf("Step=%d, time=%g, |Pi|=%g\n",i,time[1],Wave1D.l2norm(fields.u[1]))
        end
        if (mod(i,outFreq)==0)
            plt = plot(fields.grid.x, fields.u[1], xlim=xrange, ylim=yrange, label="Pi", color = :blue)
            plt = scatter!(fields.grid.x, fields.u[1], label=" ", color=:yellow)
            frame(api, plt)
            plt2 = plot(fields.grid.x, fields.u[2], xlim=xrange, ylim=(-0.1,1.0), label="Phi", color = :blue)
            plt2 = scatter!(fields.grid.x, fields.u[2], label=" ", color=:yellow)
            frame(aphi, plt2)
        end
    end

    gif(api, "aPi.gif")
    gif(aphi, "aPhi.gif")
end

function main()
  
    params = parse_commandline()
    nx = params["nx"]
    nt = params["nsteps"]
    cfl = params["cfl"]
    outFreq = params["outfreq"]

    bbox = [-10.0, 10.0]
    println("main:  nx = ",nx)
    grid = Wave1D.Grid(nx, bbox, cfl)
    println("main grid:  nx = ",grid.nx)

    fields = Wave1D.GridFields(2, grid)
    Wave1D.init_data!(fields)
    
    evolve!(fields, nt, outFreq)
end

main()

