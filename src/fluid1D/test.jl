include("grid.jl")

nx = 11
bbox = [0.0,1.0]
cfl = 0.4

g1 = Grid( nx, bbox, cfl)

f1 = GridFields(4, g1, 2)
