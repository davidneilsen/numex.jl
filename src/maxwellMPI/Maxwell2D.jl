
module Maxwell2D
using Printf

export DIM, Grid, GridFields, l2norm, rk2_step!

#c::Float64 = 1.0

const DIM=2

include("grid.jl")

include("utils.jl")

include("derivs.jl")

include("compactfd.jl")

include("rk2.jl")

include("rk4.jl")

include("eqs.jl")

include("readpars.jl")

function oops(msg::AbstractString)
    println("\n*** Oops, $msg\n")
    # close database, close files, .....
    ccall(:jl_exit, Cvoid, (Int32,), 86) # Exit code 86 to the O/S
end

end
