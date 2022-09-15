using ArgParse

include("NBody.jl")
using .NBody

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--dim", "-d"
            help = "dimension (2 or 3)1"
            arg_type = Int64
            default = 2
        "--semimajor", "-a"
            help = "binary semi-major axis"
            arg_type = Float64
            default = 2.0
        "--ecc", "-e"
            help = "binary eccentricity"
            arg_type = Float64
            default = 0.0
        "--m1"
            help = "mass 1"
            arg_type = Float64
            default = 1.0
        "--m2"
            help = "mass 2"
            arg_type = Float64
            default = 1.0
        "--m3"
            help = "mass 3"
            arg_type = Float64
            default = 1.0
        "--nb", "-n"
            help = "number of bodies"
            arg_type = Int64
            default = 2
        "time"
            help = "integration time"
            required = true
            arg_type = Float64
    end
    return parse_args(s)
end

function main()
    G = 1.0
    comline = parse_commandline()
    tf = comline["time"]
    a  = comline["semimajor"]
    e  = comline["ecc"]
    m1 = comline["m1"]
    m2 = comline["m2"]
    m3 = comline["m3"]
    nb = comline["nb"]
    ndim::Int64 = comline["dim"]

    if nb == 2
        mass = Float64[m1, m2]
    elseif nb == 3
        mass = Float64[m1, m2, m3]
    end
    pars = Params(tf, 2.0, 0.0, 4.0, 0.5)

    qb, pb = init_data_binary(ndim, m1, m2, a, e)
    bds = Bodies(qb, pb, mass, G)

    tspan = (0.0, tf)
    dt = tf/100.0
    sol = intNbdyEqs(tspan, dt, bds)
    
end

main()


