
module Maxwell2D
using Printf

export Grid, GridFields, l2norm, rk2_step!

#c::Float64 = 1.0
include("grid.jl")

include("derivs.jl")

include("rk2.jl")

include("eqs.jl")

end
