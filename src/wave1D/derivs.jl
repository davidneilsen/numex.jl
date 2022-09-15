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

