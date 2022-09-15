function diff21_x!(dxu, u, dx)
    nx, ny, nz = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    for k = 1:nz
        for j = 1:ny
            dxu[1,j,k] = (-u[1,j,k] + u[2,j,k])*idx
            for i = 2:nx-1
                dxu[i,j,k] = (u[i+1,j,k] - u[i-1,j,k]) * idx_by_2
            end
            dxu[nx,j,k] = (u[nx,j,k] - u[nx-1,j,k])*idx
        end
    end
end

function diff21_y!(dyu, u, dy)
    nx, ny, nz = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy

    for k = 1:nz
        for i = 1:nx
            dyu[i,1,k] = (-u[i,1,k] + u[i,2,k])*idy
            for j = 2:ny-1
                dyu[i,j,k] = (u[i,j+1,k] - u[i,j-1,k]) * idy_by_2
            end
            dyu[i,ny] = (u[i,ny,k] - u[i,ny-1,k])*idy
        end
    end
end

function diff21_z!(dzu, u, dz)
    nx, ny, nz = size(u)
    idz = 1.0/dz
    idz_by_2 = 0.5*idz

    for j = 1:ny
        for i = 1:nx
            dzu[i,j, 1] = (-u[i,j,1] + u[i,j,2])*idz
            for k = 2:nz-1
                dzu[i,j,k] = (u[i,j,k+1,k] - u[i,j,k-1]) * idz_by_2
            end
            dzu[i,j,nz] = (u[i,j,nz] - u[i,j,nz-1])*idz
        end
    end
end



function diff22_x!(dxu, u, dx)
    nx, ny, nz = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    for k = 1:nz
        for j = 1:ny
            dxu[1,j,k] = (-3.0*u[1,j,k] + 4.0*u[2,j,k] - u[3,j,k]) * idx_by_2
            for i = 2:nx-1
                dxu[i,j,k] = (u[i+1,j,k] - u[i-1,j,k]) * idx_by_2
            end
            dxu[nx,j,k] = (u[nx-2,j,k] - 4.0*u[nx-1,j,k] + 3.0*u[nx,j,k]) * idx_by_2
        end
    end
end

function diff22_y!(dyu, u, dy)
    nx, ny, nz = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy

    for k = 1:nz
        for i = 1:nx
            dyu[i,1,k] = (-3.0*u[i,1,k] + 4.0*u[i,2,k] - u[i,3,k]) * idy_by_2
            for j = 2:ny-1
                dyu[i,j,k] = (u[i,j+1,k] - u[i,j-1,k]) * idy_by_2
            end
            dyu[i,ny,k] = (u[i,ny-2,k] - 4.0*u[i,ny-1,k] + 3.0*u[i,ny,k]) * idy_by_2
        end
    end
end

function diff22_z!(dzu, u, dz)
    nx, ny, nz = size(u)
    idz = 1.0/dz
    idz_by_2 = 0.5*idz

    for j = 1:ny
        for i = 1:nx
            dzu[i,j,1] = (-3.0*u[i,j,1] + 4.0*u[i,j,2] - u[i,j,3]) * idz_by_2
            for k = 2:nz-1
                dzu[i,j,k] = (u[i,j,k+1] - u[i,j,k-1]) * idz_by_2
            end
            dzu[i,j,nz] = (u[i,j,nz-2] - 4.0*u[i,j,nz-1] + 3.0*u[i,j,nz]) * idz_by_2
        end
    end
end


function diff42_x!(dxu, u, dx)
    nx, ny, nz = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx
    idx_by_12 = idx/12.0

    for k = 1:nz
        for j = 1:ny
            dxu[1,j,k] = (-3.0*u[1,j,k] + 4.0*u[2,j,k] - u[3,j,k]) * idx_by_2
            dxu[2,j,k] = (-u[1,j,k] + u[3,j,k]) * idx_by_2
            for i = 3:nx-2
                dxu[i,j,k] = (u[i-2,j,k] - 8.0*u[i-1,j,k] + 8.0*u[i+1,j,k] - u[i+2,j,k])*idx_by_12
            end
            dxu[nx-1,j,k] = ( - u[nx-2,j,k] + u[nx,j,k]) * idx_by_2
            dxu[nx,j,k] = (u[nx-2,j,k] - 4.0*u[nx-1,j,k] + 3.0*u[nx,j,k]) * idx_by_2
        end
    end
end

function diff42_y!(dyu, u, dy)
    nx, ny, nz = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy
    idy_by_12 = idy/12.0

    for k = 1:nz
        for i = 1:nx
            dyu[i,1,k] = (-3.0*u[i,1,k] + 4.0*u[i,2,k] - u[i,3,k]) * idy_by_2
            dyu[i,2,k] = (-u[i,1,k] + u[i,3,k]) * idy_by_2
            for j = 3:ny-2
                dyu[i,j,k] = (u[i,j-2,k] - 8.0*u[i,j-1,k] + 8.0*u[i,j+1,k] - u[i,j+2,k])*idy_by_12
            end
            dyu[i,ny-1,k] = ( - u[i,ny-2,k] + u[i,ny,k]) * idy_by_2
            dyu[i,ny,k] = (u[i,ny-2,k] - 4.0*u[i,ny-1,k] + 3.0*u[i,ny,k]) * idy_by_2
        end
    end
end

function diff42_z!(dzu, u, dz)
    nx, ny, nz = size(u)
    idz = 1.0/dy
    idz_by_2 = 0.5*idz
    idz_by_12 = idz/12.0

    for j = 1:ny
        for i = 1:nx
            dyu[i,j,1] = (-3.0*u[i,j,1] + 4.0*u[i,j,2] - u[i,j,3]) * idz_by_2
            dyu[i,j,2] = (-u[i,j,1] + u[i,j,3]) * idz_by_2
            for k = 3:nz-2
                dyu[i,j,j] = (u[i,j,j-2] - 8.0*u[i,j,j-1] + 8.0*u[i,j,j+1] - u[i,j,j+2])*idz_by_12
            end
            dyu[i,j,nz-1] = ( - u[i,j,nz-2] + u[i,j,nz]) * idz_by_2
            dyu[i,j,nz] = (u[i,j,nz-2] - 4.0*u[i,j,nz-1] + 3.0*u[i,j,nz]) * idz_by_2
        end
    end
end


function kodiss!(f, u, dx, dy, dz)
    kodiss_x!(f, u, dx)
    kodiss_y!(f, u, dy)
    kodiss_z!(f, u, dz)
end

function kodiss_x!(f, u, dx)
    nx, ny, nz = size(f)
    sigma = 0.05
    idx = sigma/dx
    for k = 1:nz
        for j = 1:ny
            f[1,j,k] -= 2.0*(u[3,j,k] - 2.0*u[2,j,k] + u[1,j,k])*idx
            f[2,j,k] -= (u[4,j,k] - 4.0*u[3,j,k] + 5.0*u[2,j,k] - 2.0*u[1,j,k])*idx
            for i = 3:nx-2
                f[i,j,k] -= (u[i+2,j,k] - 4.0*u[i+1,j,k] + 6.0*u[i,j,k]
                                    - 4.0*u[i-1,j,k] + u[i-2,j,k])*idx
            end                     
            f[nx-1,j,k] -= (u[nx-3,j,k] - 4.0*u[nx-2,j,k] + 5.0*u[nx-1,j,k] - 2.0*u[nx,j,k]  ) *idx
            f[nx,j,k] -= 2.0 * (u[nx-2,j,k] - 2.0*u[nx-1,j,k] + u[nx,j,k])*idx
        end
    end

end

function kodiss_y!(f, u, dy)
    nx, ny, nz = size(f)
    sigma = 0.05
    idy = sigma/dy

    for k = 1:nz
        for i = 1:nx
            f[i,1,k] -= 2.0*(u[i,3,k] - 2.0*u[i,2,k] + u[i,1,k])*idy
            f[i,2,k] -= (u[i,4,k] - 4.0*u[i,3,k] + 5.0*u[i,2,k] - 2.0*u[i,1,k]  )*idy
            for j = 3:ny-2
                f[i,j,k] -= (u[i,j+2,k] - 4.0*u[i,j+1,k] + 6.0*u[i,j,k]
                                    - 4.0*u[i,j-1,k] + u[i,j-2,k])*idy
            end
            f[i,ny-1,k] -= (u[i,ny-3,k] - 4.0*u[i,ny-2,k] + 5.0*u[i,ny-1,k] - 2.0*u[i,ny,k])*idy
            f[i,ny,k] -= 2.0*(u[i,ny-2,k] - 2.0*u[i,ny-1,k] + u[i,ny,k]  )*idy
        end
    end

end

function kodiss_z!(f, u, dz)
    nx, ny, nz = size(f)
    sigma = 0.05
    idz = sigma/dz

    for j = 1:ny
        for i = 1:nx
            f[i,j,1] -= 2.0*(u[i,j,3] - 2.0*u[i,j,2] + u[i,j,1])*idz
            f[i,j,2] -= (u[i,j,4] - 4.0*u[i,j,3] + 5.0*u[i,j,2] - 2.0*u[i,j,1]  )*idz
            for k = 3:nz-2
                f[i,j,k] -= (u[i,j,k+2] - 4.0*u[i,j,k+1] + 6.0*u[i,j,k]
                                    - 4.0*u[i,j,k-1] + u[i,j,k-2])*idz
            end
            f[i,j,nz-1] -= (u[i,j,nz-3] - 4.0*u[i,j,nz-2] + 5.0*u[i,j,nz-1] - 2.0*u[i,j,nz])*idz
            f[i,j,nz] -= 2.0*(u[i,j,nz-2] - 2.0*u[i,j,nz-1] + u[i,j,nz]  )*idz
        end
    end

end

