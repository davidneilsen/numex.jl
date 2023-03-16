using Printf
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nx", "-n"
            help = "number of points in each direction"
            arg_type = Int
            default = 51
        "--cfl", "-c"
            help = "Courant number"
            arg_type = Float64
            default = 0.1
        "--idtype"
            help = "Initial Data type"
            arg_type = Int
            default = 0
        "nsteps"
            help = "Number of steps"
            arg_type = Int
            required = true
    end
    return parse_ags(s)
end

function main()
    params = parse_commandline()
    set_initial!(fields)

end

main()
