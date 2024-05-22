
#=---------------------------------------------------------------------
 =
 =  This routine sets the initial data for the wave equation.
 =  Currently it sets time-symmetric initial data, where the time
 =  derivative of the scalar field is set to zero.  The spatial 
 =  derivative is set to a Gaussian.
 =
 =   Pi  = \partial_t \phi
 =   Phi = \partial_x \phi
 =
 =---------------------------------------------------------------------=#
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

#=---------------------------------------------------------------------
 =
 =  The evolution equations for the scalar field.  This calculates the RHS
 =  for the Method of Lines.
 =
 =---------------------------------------------------------------------=#
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

    # calculate spatial derivatives of the Pi and Phi
    diff22_x!(dxPi, Pi, dx)
    diff22_x!(dxPhi, Phi, dx)

    # out-going Sommerfeld boundary condition on the left
    dtPi[1] = idx*(Pi[2] - Pi[1])
    dtPhi[1] = idx*(Phi[2] - Phi[1])

    # evaluate the RHS of differential equations in the grid interior
    for i = 2:nx-1
        dtPi[i] = dxPhi[i]
        dtPhi[i] = dxPi[i]
    end

    # out-going Sommerfeld boundary condition on the right
    dtPi[nx] = -idx*(Pi[nx] - Pi[nx-1])
    dtPhi[nx] = -idx*(Phi[nx] - Phi[nx-1])

end

