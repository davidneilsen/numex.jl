module CompactFDWave1D
using Printf
using LinearAlgebra

export Grid, GridFields, gauss_init_data!, l2norm

include("grid.jl")

include("compactderivs.jl")

include("rk4.jl")

include("eqs.jl")

include("derivs.jl")

end
