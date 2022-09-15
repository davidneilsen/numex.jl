
function init_data!(fields)
    nx, ny, nz = fields.grid.shp

    Ex = fields.u[1]
    Ey = fields.u[2]
    x = fields.grid.x
    y = fields.grid.y
    z = fields.grid.z

    amp1 = 1.0
    lambda1 = 1.0

    for k=1:nz, j=1:ny, i=1:nx
        r = sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        Ephi = - 8.0*amp1*lambda1*lambda1*exp(-lambda1*r*r)
        Ex[i,j,k] = - y[j] * Ephi
        Ey[i,j,k] =   x[i] * Ephi
    end

end

function maxwellEqs!(dtu, u, dxu, dyu, dzu, x, y, z, dx, dy, dz, time)
    dtEx = dtu[1]
    dtEy = dtu[2]
    dtEz = dtu[3]
    dtBx = dtu[4]
    dtBy = dtu[5]
    dtBz = dtu[6]

    dxEy = dxu[2]
    dxEz = dxu[3]
    dxBy = dxu[5]
    dxBz = dxu[6]

    dyEx = dyu[1]
    dyEz = dyu[3]
    dyBx = dyu[4]
    dyBz = dyu[6]

    dzEx = dzu[1]
    dzEy = dzu[2]
    dzBx = dzu[4]
    dzBy = dzu[5]

    @. diff22_x!(dxu, u, dx)
    @. diff22_y!(dyu, u, dy)
    @. diff22_z!(dzu, u, dz)

    @. dtEx = dyBz - dzBy
    @. dtEy = dzBx - dxBz
    @. dtEy = dxBy - dyBx

    @. dtBx = dzEy - dyEz
    @. dtBy = dxEz - dzEx
    @. dtBz = dyEx - dxEy

    for i = 1:length(u)
        sommerfeld_bcs(dtu[i], u[i], dxu[i], dyu[i], dzu[i], x, y, z)
    end

end

function sommerfeld_bcs(dtu, u, dxu, dyu, dzu, x, y, z)

    nx, ny, nz = size(u)

    u0::Float64 = 0.0
    u_falloff::Float64 = 2.0

    ############  j = 1
    j = 1
    for k = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  j = ny
    j = ny
    for k = 1:nz
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  i = 1
    i = 1
    for k = 1:nz
        for j = 1:ny
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  i = nx
    i = nx
    for k = 1:nz
        for j = 1:ny
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  k = 1
    k = 1
    for j = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  k = nz
    k = nz
    for j = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] + z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end
    
end
