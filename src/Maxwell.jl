
module Maxwell
using Printf

export Grid, GridFields, init_data!, l2norm

#c::Float64 = 1.0
include("maxwell/grid.jl")

include("maxwell/derivs.jl")

include("maxwell/rk2.jl")

include("maxwell/eqs.jl")

end
