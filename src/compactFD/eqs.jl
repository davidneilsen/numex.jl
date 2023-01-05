
function gauss_init_data!(f, x)
    nx = length(x)

    amp = 1.0
    lambda = 1.0
    x0 = 0.0

    @. f = amp*exp(-((x - x0)/lambda )^2)

end

function waveEqs!(dtu, u, dxu, x, dx, time, A, B, order, method)
    nx = length(x)
    idx = 1.0/dx
    idx_by_2 = 1.0/(2.0*dx)

    dtPi = dtu[1]
    dtPhi = dtu[2]
    Pi = u[1]
    Phi = u[2]
    dxPi = dxu[1]
    dxPhi = dxu[2]
    
    if method == 1
        dxPi = A \ (B*Pi)
        dxPhi = A \ (B*Phi)
    else
        #print("calling FD derivs")
        if order == 4
            diff44_x!(dxPi, Pi, dx)
            diff44_x!(dxPhi, Phi, dx)
        else
            diff666_x!(dxPi, Pi, dx)
            diff666_x!(dxPhi, Phi, dx)
        end
    end

    dtPi[1] = dxPi[1]
    dtPhi[1] = dxPhi[1]

    for i = 2:nx-1
        dtPi[i] = dxPhi[i]
        dtPhi[i] = dxPi[i]
    end

    dtPi[nx] = - dxPi[nx]
    dtPhi[nx] = - dxPhi[nx]

end

function advectionEq!(dtu, u, dxu, x, dx, time, A, B)
    nx = length(x)
    idx = 1.0/dx
    idx_by_2 = 1.0/(2.0*dx)

    dtPi = dtu[1]
    Pi = u[1]
    dxPi = dxu[1]

    dxPi = A \ (B*Pi)

    for i = 1:nx-1
        dtPi[i] = dxPi[i]
    end

    dtPi[nx] = -dxPi[nx]
end



