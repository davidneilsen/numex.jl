
module Maxwell
using Printf

export Grid, GridFields, init_data!, l2norm

#c::Float64 = 1.0
include("grid.jl")

include("derivs.jl")

include("rk2.jl")

include("eqs.jl")

end
