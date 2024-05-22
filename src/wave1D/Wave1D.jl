
module Wave1D
using Printf

export Grid, GridFields, init_data!, l2norm, output_curve

include("grid.jl")

include("derivs.jl")

include("rk2.jl")

include("eqs.jl")

include("curve.jl")

end
