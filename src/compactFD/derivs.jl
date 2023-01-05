function diff21_x!(dxu, u, dx)
    nx = length(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    dxu[1] = (-u[1] + u[2])*idx
    for i = 2:nx-1
        dxu[i] = (u[i+1] - u[i-1]) * idx_by_2
    end
    dxu[nx] = (u[nx] - u[nx-1])*idx
end

function diff22_x!(dxu, u, dx)
    nx = length(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    dxu[1] = (-3.0*u[1] + 4.0*u[2] - u[3]) * idx_by_2
    for i = 2:nx-1
        dxu[i] = (u[i+1] - u[i-1]) * idx_by_2
    end
    dxu[nx] = (u[nx-2] - 4.0*u[nx-1] + 3.0*u[nx]) * idx_by_2
end

function diff42_x!(dxu, u, dx)
    nx = length(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx
    idx_by_12 = idx/12.0

    dxu[1] = (-3.0*u[1] + 4.0*u[2] - u[3]) * idx_by_2
    dxu[2] = (-u[1] + u[3]) * idx_by_2
    for i = 3:nx-2
        dxu[i] = (u[i-2] - 8.0*u[i-1] + 8.0*u[i+1] - u[i+2])*idx_by_12
    end
    dxu[nx-1] = ( - u[nx-2] + u[nx]) * idx_by_2
    dxu[nx] = (u[nx-2] - 4.0*u[nx-1] + 3.0*u[nx]) * idx_by_2
end

function diff44_x!(dxu, u, dx)
    nx = length(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx
    idx_by_12 = idx/12.0

    dxu[1] = ( - 25.0 * u[1  ] + 48.0 * u[2] - 36.0 * u[3] + 16.0 * u[4] -  3.0 * u[5] ) * idx_by_12
    dxu[2] = (  -  3.0 * u[1] - 10.0 * u[2] + 18.0 * u[3] -  6.0 * u[4] + u[5] ) * idx_by_12
    for i = 3:nx-2
        dxu[i] = (u[i-2] - 8.0*u[i-1] + 8.0*u[i+1] - u[i+2])*idx_by_12
    end

    i=nx-1
    dxu[i] = (  -  u[i-3] +  6.0 * u[i-2] - 18.0 * u[i-1] + 10.0 * u[i  ] +  3.0 * u[i+1] ) * idx_by_12
    i=nx
    dxu[i] = ( 3.0 * u[i-4] - 16.0 * u[i-3] + 36.0 * u[i-2] - 48.0 * u[i-1] + 25.0 * u[i] ) * idx_by_12

end

function diff666_x!(dxu, u, dx)
    nx  = length(u)
    idx = 1.0/dx
    idx_by_60 = idx / 60.0
    
    dxu[1] = ( - 147.0*u[1] + 360.0*u[2] - 450.0*u[3] + 400.0*u[4] - 225.0*u[5] +  72.0*u[6] -  10.0*u[7] ) * idx_by_60
    dxu[2] = ( -  10.0 * u[1] - 77.0 * u[2] + 150.0 * u[3] - 100.0*u[4] + 50.0*u[5] - 15.0*u[6] + 2.0*u[7] ) * idx_by_60
    dxu[3] = ( 2.0*u[1] - 24.0 *u[2] - 35.0*u[3] + 80.0*u[4] - 30.0*u[5] + 8.0*u[6] - u[7] ) * idx_by_60

    for i = 4:nx-3
        dxu[i] = ( -u[i-3] + 9.0*u[i-2] - 45.0*u[i-1] + 45.0*u[i+1] - 9.0*u[i+2] + u[i+3] ) * idx_by_60
    end

    i = nx-2
    dxu[i] = ( u[i-4] - 8.0 * u[i-3] + 30.0 * u[i-2] - 80.0 * u[i-1] + 35.0 * u[i  ] + 24.0 * u[i+1] -  2.0 * u[i+2] ) * idx_by_60
    i=nx-1
    dxu[i] = (  -   2.0 * u[i-5] + 15.0 * u[i-4] - 50.0 * u[i-3] + 100.0 * u[i-2] - 150.0 * u[i-1] +  77.0 * u[i  ] +  10.0 * u[i+1] ) * idx_by_60
    i=nx
    dxu[i] = ( 10.0 * u[i-6] -  72.0 * u[i-5] + 225.0 * u[i-4] - 400.0 * u[i-3] + 450.0 * u[i-2] - 360.0 * u[i-1] + 147.0 * u[i  ] ) * idx_by_60

end

function kodiss!(f, u, dx)

    nx = length(f)
    sigma = 0.05
    idx = sigma/dx

    f[1] -= 2.0*(u[3] - 2.0*u[2] + u[1])*idx
    f[2] -= (u[4] - 4.0*u[3] + 5.0*u[2] - 2.0*u[1])*idx
    for i = 3:nx-2
        f[i] -= (u[i+2] - 4.0*u[i+1] + 6.0*u[i]
                            - 4.0*u[i-1] + u[i-2])*idx
    end
    f[nx-1] -= (u[nx-3] - 4.0*u[nx-2] + 5.0*u[nx-1] - 2.0*u[nx]  ) *idx
    f[nx] -= 2.0 * (u[nx-2] - 2.0*u[nx-1] + u[nx])*idx

end

