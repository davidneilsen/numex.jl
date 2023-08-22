####################################################################
#
#  Finite difference operators.
#  Second and fourth order
#
####################################################################

function diff21_x!(dxu, u, dx)
    nx, ny = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    for j = 1:ny
        dxu[1,j] = (-u[1,j] + u[2,j])*idx
        for i = 2:nx-1
            dxu[i,j] = (u[i+1,j] - u[i-1,j]) * idx_by_2
        end
        dxu[nx,j] = (u[nx,j] - u[nx-1,j])*idx
    end
end

function diff21_y!(dyu, u, dy)
    nx, ny = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy

    for i = 1:nx
        dyu[i,1] = (-u[i,1] + u[i,2])*idy
        for j = 2:ny-1
            dyu[i,j] = (u[i,j+1] - u[i,j-1]) * idy_by_2
        end
        dyu[i,ny] = (u[i,ny] - u[i,ny-1])*idy
    end
end


function diff22_x!(dxu, u, dx)
    nx, ny = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx

    for j = 1:ny
        dxu[1,j] = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx_by_2
        for i = 2:nx-1
            dxu[i,j] = (u[i+1,j] - u[i-1,j]) * idx_by_2
        end
        dxu[nx,j] = (u[nx-2,j] - 4.0*u[nx-1,j] + 3.0*u[nx,j]) * idx_by_2
    end
end

function diff22_y!(dyu, u, dy)
    nx, ny = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy

    for i = 1:nx
        dyu[i,1] = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy_by_2
        for j = 2:ny-1
            dyu[i,j] = (u[i,j+1] - u[i,j-1]) * idy_by_2
        end
        dyu[i,ny] = (u[i,ny-2] - 4.0*u[i,ny-1] + 3.0*u[i,ny]) * idy_by_2
    end
end

function diff42_x!(dxu, u, dx)
    nx, ny = size(u)
    idx = 1.0/dx
    idx_by_2 = 0.5*idx
    idx_by_12 = idx/12.0

    for j = 1:ny
        dxu[1,j] = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx_by_2
        dxu[2,j] = (-u[1,j] + u[3,j]) * idx_by_2
        for i = 3:nx-2
            dxu[i,j] = (u[i-2,j] - 8.0*u[i-1,j] + 8.0*u[i+1,j] - u[i+2,j])*idx_by_12
        end
        dxu[nx-1,j] = ( - u[nx-2,j] + u[nx,j]) * idx_by_2
        dxu[nx,j] = (u[nx-2,j] - 4.0*u[nx-1,j] + 3.0*u[nx,j]) * idx_by_2
    end
end

function diff42_y!(dyu, u, dy)
    nx, ny = size(u)
    idy = 1.0/dy
    idy_by_2 = 0.5*idy
    idy_by_12 = idy/12.0

    for i = 1:nx
        dyu[i,1] = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy_by_2
        dyu[i,2] = (-u[i,1] + u[i,3]) * idy_by_2
        for j = 3:ny-2
            dyu[i,j] = (u[i,j-2] - 8.0*u[i,j-1] + 8.0*u[i,j+1] - u[i,j+2])*idy_by_12
        end
        dyu[i,ny-1] = ( - u[i,ny-2] + u[i,ny]) * idy_by_2
        dyu[i,ny] = (u[i,ny-2] - 4.0*u[i,ny-1] + 3.0*u[i,ny]) * idy_by_2
    end
end

function diff666_x!(Dxu, u, dx)
    nx, ny = size(u)
    idx = 1.00/dx
    idx_by_60 = idx / 60.0
    for j = 1:ny
        i=1
        Dxu[i,j] = ( - 147.0 * u[i,j]     
                        + 360.0 * u[i+1,j]     
                        - 450.0 * u[i+2,j]     
                        + 400.0 * u[i+3,j]     
                        - 225.0 * u[i+4,j]     
                        +  72.0 * u[i+5,j]     
                        -  10.0 * u[i+6,j]     
                      ) * idx_by_60
        i=2
        Dxu[2,j] = (  -  10.0 * u[i-1,j]     
                         -  77.0 * u[i  ,j]     
                         + 150.0 * u[i+1,j]     
                         - 100.0 * u[i+2,j]     
                         +  50.0 * u[i+3,j]     
                         -  15.0 * u[i+4,j]     
                         +   2.0 * u[i+5,j]     
                      ) * idx_by_60
        i=3
        Dxu[3,j] = (      2.0 * u[i-2,j]     
                         -  24.0 * u[i-1,j]     
                         -  35.0 * u[i  ,j]     
                         +  80.0 * u[i+1,j]     
                         -  30.0 * u[i+2,j]     
                         +   8.0 * u[i+3,j]     
                         -          u[i+4,j]     
                      ) * idx_by_60
        for i = 4:nx-3
            Dxu[i,j] = ( -         u[i-3,j]     
                          +  9.0 * u[i-2,j]     
                          - 45.0 * u[i-1,j]     
                          + 45.0 * u[i+1,j]     
                          -  9.0 * u[i+2,j]     
                          +         u[i+3,j]     
                        ) * idx_by_60
        end
        i=nx-2
        Dxu[i,j] = (            u[i-4,j]     
                         -  8.0 * u[i-3,j]     
                         + 30.0 * u[i-2,j]     
                         - 80.0 * u[i-1,j]     
                         + 35.0 * u[i  ,j]     
                         + 24.0 * u[i+1,j]     
                         -  2.0 * u[i+2,j]     
                      ) * idx_by_60
        i=nx-1
        Dxu[i,j] = (  -   2.0 * u[i-5,j]     
                         +  15.0 * u[i-4,j]     
                         -  50.0 * u[i-3,j]     
                         + 100.0 * u[i-2,j]     
                         - 150.0 * u[i-1,j]     
                         +  77.0 * u[i  ,j]     
                         +  10.0 * u[i+1,j]     
                       ) * idx_by_60
        i=nx
        Dxu[i,j] = (     10.0 * u[i-6,j]     
                         -  72.0 * u[i-5,j]     
                         + 225.0 * u[i-4,j]     
                         - 400.0 * u[i-3,j]     
                         + 450.0 * u[i-2,j]     
                         - 360.0 * u[i-1,j]     
                         + 147.0 * u[i,j]     
                       ) * idx_by_60
    end

end

function diff666_y!(Dyu, u, dy)
    nx, ny = size(u)
    idy = 1.00/dy
    idy_by_60 = idy / 60.0
    for i = 1:nx
        j=1
        Dyu[i,j] = ( - 147.0 * u[i,j]     
                        + 360.0 * u[i,j+1]     
                        - 450.0 * u[i,j+2]     
                        + 400.0 * u[i,j+3]     
                        - 225.0 * u[i,j+4]     
                        +  72.0 * u[i,j+5]     
                        -  10.0 * u[i,j+6]     
                      ) * idy_by_60
        j=2
        Dyu[2] = (  -  10.0 * u[i,j-1]     
                         -  77.0 * u[i,j  ]     
                         + 150.0 * u[i,j+1]     
                         - 100.0 * u[i,j+2]     
                         +  50.0 * u[i,j+3]     
                         -  15.0 * u[i,j+4]     
                         +   2.0 * u[i,j+5]     
                      ) * idy_by_60
        j=3
        Dyu[3] = (      2.0 * u[i,j-2]     
                         -  24.0 * u[i,j-1]     
                         -  35.0 * u[i,j  ]     
                         +  80.0 * u[i,j+1]     
                         -  30.0 * u[i,j+2]     
                         +   8.0 * u[i,j+3]     
                         -          u[i,j+4]     
                      ) * idy_by_60
        for j = 4:ny-3
            Dyu[i,j] = ( -         u[i,j-3]     
                          +  9.0 * u[i,j-2]     
                          - 45.0 * u[i,j-1]     
                          + 45.0 * u[i,j+1]     
                          -  9.0 * u[i,j+2]     
                          +        u[i,j+3]     
                        ) * idy_by_60
        end
        j=ny-2
        Dyu[i,j] = (            u[i,j-4]     
                         -  8.0 * u[i,j-3]     
                         + 30.0 * u[i,j-2]     
                         - 80.0 * u[i,j-1]     
                         + 35.0 * u[i,j  ]     
                         + 24.0 * u[i,j+1]     
                         -  2.0 * u[i,j+2]     
                      ) * idy_by_60
        j=ny-1
        Dyu[i,j] = (  -   2.0 * u[i,j-5]     
                         +  15.0 * u[i,j-4]     
                         -  50.0 * u[i,j-3]     
                         + 100.0 * u[i,j-2]     
                         - 150.0 * u[i,j-1]     
                         +  77.0 * u[i,j  ]     
                         +  10.0 * u[i,j+1]     
                       ) * idy_by_60
        j=ny
        Dyu[i,j] = (     10.0 * u[i,j-6]     
                         -  72.0 * u[i,j-5]     
                         + 225.0 * u[i,j-4]     
                         - 400.0 * u[i,j-3]     
                         + 450.0 * u[i,j-2]     
                         - 360.0 * u[i,j-1]     
                         + 147.0 * u[i,j]     
                       ) * idy_by_60
    end

end


function kodiss_nobndry!(f, u, dx, dy)
    nx, ny = size(f)
    sigma = 0.05
    idx = sigma/dx
    idy = sigma/dy
    for j = 3:ny-2, i = 3:nx-2
        f[i,j] -= (u[i+2,j] - 4.0*u[i+1,j] + 6.0*u[i,j]
                            - 4.0*u[i-1,j] + u[i-2,j])*idx
                + (u[i,j+2] - 4.0*u[i,j+1] + 6.0*u[i,j]
                            - 4.0*u[i,j-1] + u[i,j-2])*idy
    end
end

function kodiss!(f, u, dx, dy)
    kodiss_x!(f, u, dx)
    kodiss_y!(f, u, dy)
end

function kodiss_x!(f, u, dx)
    nx, ny = size(f)
    sigma = 0.05
    idx = sigma/dx
    for j = 1:ny
        f[1,j] -= 2.0*(u[3,j] - 2.0*u[2,j] + u[1,j])*idx
        f[2,j] -= (u[4,j] - 4.0*u[3,j] + 5.0*u[2,j] - 2.0*u[1,j])*idx
        for i = 3:nx-2
            f[i,j] -= (u[i+2,j] - 4.0*u[i+1,j] + 6.0*u[i,j]
                                - 4.0*u[i-1,j] + u[i-2,j])*idx
        end
        f[nx-1,j] -= (u[nx-3,j] - 4.0*u[nx-2,j] + 5.0*u[nx-1,j] - 2.0*u[nx,j]  ) *idx
        f[nx,j] -= 2.0 * (u[nx-2,j] - 2.0*u[nx-1,j] + u[nx,j])*idx
    end

end

function kodiss_y!(f, u, dy)
    nx, ny = size(f)
    sigma = 0.05
    idy = sigma/dy

    for i = 1:nx
        f[i,1] -= 2.0*(u[i,3] - 2.0*u[i,2] + u[i,1])*idy
        f[i,2] -= (u[i,4] - 4.0*u[i,3] + 5.0*u[i,2] - 2.0*u[i,1]  )*idy
        for j = 3:ny-2
            f[i,j] -= (u[i,j+2] - 4.0*u[i,j+1] + 6.0*u[i,j]
                                - 4.0*u[i,j-1] + u[i,j-2])*idy
        end
        f[i,ny-1] -= (u[i,ny-3] - 4.0*u[i,ny-2] + 5.0*u[i,ny-1] - 2.0*u[i,ny])*idy
        f[i,ny] -= 2.0*(u[i,ny-2] - 2.0*u[i,ny-1] + u[i,ny]  )*idy
    end
 
end

function cdiff_x!(dxu, u, Dx, u1d, dxu1d)
    nx, ny = size(u)
    for j = 1:ny
        for i = 1:nx
            u1d[i] = u[i,j]
        end
        #mul!(dxu1d, Dx, u1d)
        dxu1d .= 0.0
        for ii = 1:nx, jj = 1:nx
            dxu1d[ii] += Dx[ii,jj]*u1d[jj]
        end 
        for i = 1:nx
            dxu[i,j] = dxu1d[i]
        end
    end
end

function cdiff_y!(dyu, u, Dy, u1d, dyu1d)
    nx, ny = size(u)
    for i = 1:nx
        for j = 1:ny
            u1d[j] = u[i,j]
        end
        #mul!(dyu1d, Dy, u1d)
        dyu1d .= 0.0
        for ii = 1:ny
            for jj = 1:ny
                dyu1d[ii] += Dy[ii,jj]*u1d[jj]
            end
        end
        for j = 1:ny
            dyu[i,j] = dyu1d[j]
        end
    end
end

function cfilter!(u, dvars)
    #println("...applying cfilter!")
    cfilter_x!(u, dvars.Fx1d, dvars.u1x, dvars.du1x)
    cfilter_y!(u, dvars.Fy1d, dvars.u1y, dvars.du1y)
end

function cfilter_x!(u, F, uf1d, u1d)
    nx, ny = size(u)
    for j = 1:ny
        for i = 1:nx
            u1d[i] = u[i,j]
        end
        #mul!(dxu1d, Dx, u1d)
        for ii = 1:nx
            uf1d[ii] = 0.0
            for jj = 1:nx
                uf1d[ii] += F[ii,jj]*u1d[jj]
            end 
        end 
        for i = 1:nx
            u[i,j] += uf1d[i]
        end
    end
end

function cfilter_y!(u, F, uf1d, u1d)
    nx, ny = size(u)
    for i = 1:nx
        for j = 1:ny
            u1d[j] = u[i,j]
        end
        #mul!(dyu1d, Dy, u1d)
        for ii = 1:ny
            uf1d[ii] = 0.0
            for jj = 1:ny
                uf1d[ii] += F[ii,jj]*u1d[jj]
            end
        end
        for j = 1:ny
            u[i,j] += uf1d[j]
        end
    end
end



