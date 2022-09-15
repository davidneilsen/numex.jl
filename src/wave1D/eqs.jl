
function init_data!(fields)
    nx = fields.grid.nx

    Pi = fields.u[1]
    Phi = fields.u[2]
    x = fields.grid.x

    amp = 1.0
    lambda = 1.0
    x0 = 0.0

    for i=1:nx
        Pi[i] = 0.0
        Phi[i] = amp*exp(-( (x[i]-x0)/lambda )^2)
    end

end

function waveEqs!(dtu, u, dxu, x, dx, time)
    nx = length(x)
    idx = 1.0/dx
    idx_by_2 = 1.0/(2.0*dx)

    dtPi = dtu[1]
    dtPhi = dtu[2]
    Pi = u[1]
    Phi = u[2]
    dxPi = dxu[1]
    dxPhi = dxu[2]

    diff22_x!(dxPi, Pi, dx)
    diff22_x!(dxPhi, Phi, dx)

    dtPi[1] = idx*(Pi[2] - Pi[1])
    dtPhi[1] = idx*(Pi[2] - Pi[1])

    for i = 2:nx-1
        dtPi[i] = dxPhi[i]
        dtPhi[i] = dxPi[i]
    end

    dtPi[nx] = -idx*(Pi[nx] - Pi[nx-1])
    dtPhi[nx] = -idx*(Phi[nx] - Phi[nx-1])

end

