
module Maxwell2D
using Printf

export DIM, Grid, GridFields, l2norm, rk2_step!

#c::Float64 = 1.0

const DIM=2

include("grid.jl")

include("utils.jl")

include("derivs.jl")

include("rk2.jl")

include("eqs.jl")

include("readpars.jl")


end
