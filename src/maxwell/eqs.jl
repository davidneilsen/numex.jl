
function init_data!(fields)
    nx = fields.grid.nx
    ny = fields.grid.ny

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]
    x = fields.grid.x
    y = fields.grid.y

    amp1 = 1.0
    lambda1 = 1.0

    for j=1:ny, i=1:nx
        r = sqrt(x[i]*x[i] + y[j]*y[j])
        Bz[i,j] = - 8.0*amp1*lambda1*lambda1*exp(-lambda1*r*r)
    end

end

function Xinit_data!(fields)
    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    md = 2.0
    nd = 2.0
    c = 1.0

    dimx = fields.grid.xmax - fields.grid.xmin
    dimy = fields.grid.ymax - fields.grid.ymin
    nx = fields.grid.nx
    ny = fields.grid.ny
    dx = fields.grid.dx
    dy = fields.grid.dy
    dt = fields.grid.dt
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]

    println("init_data:  nx = ",nx," ny = ",ny)
    for j=1:ny, i=1:nx
        Bz[i,j] = (- cos(md*pi*((i-0.5)*dx/dimx)) 
                        * cos(nd*pi*((j-0.5)*dy/dimy))
                        * cos(omega*(-0.5*dt)) )
#        @printf("i=%d, j=%d, Ex=%g, Ey=%g, Bz = %g",i,j,Ex[i,j],Ey[i,j],Bz[i,j])
    end  
end

function maxwell_TE!(dtu, u, dxu, dyu, x, y, dx, dy, time)
    dtEx = dtu[1]
    dtEy = dtu[2]
    dtHz = dtu[3]
    Ex = u[1]
    Ey = u[2]
    Hz = u[3]
    dyEx = dyu[1]
    dxEy = dxu[2]
    dxHz = dxu[3]
    dyHz = dyu[3]

    diff22_y!(dyHz, Hz, dy)
    diff22_x!(dxHz, Hz, dx)
    diff22_x!(dxEy, Ey, dy)
    diff22_y!(dyEx, Ex, dx)

    @. dtEx = dyHz
    @. dtEy = -dxHz
    @. dtHz = dyEx - dxEy

#
    sommerfeld_bcs(dtu[1], u[1], x, y)
    sommerfeld_bcs(dtu[2], u[2], x, y)
    sommerfeld_bcs(dtu[3], u[3], x, y)
#=
    sommerfeld_bcs_deriv(dtu[1], u[1], dxu[1], dyu[1], x, y)
    sommerfeld_bcs_deriv(dtu[2], u[2], dxu[2], dyu[2], x, y)
    sommerfeld_bcs_deriv(dtu[3], u[3], dxu[3], dyu[3], x, y)
=#
#=
    dirichlet_bcs(dtu[1], u[1], x, y)
    dirichlet_bcs(dtu[2], u[2], x, y)
    dirichlet_bcs(dtu[3], u[3], x, y)
=#

end

function maxwell_TM!(dtu, u, dxu, dyu, x, y, dx, dy, time)
    dtHx = dtu[1]
    dtHy = dtu[2]
    dtEz = dtu[3]
    Hx = u[1]
    Hy = u[2]
    Ez = u[3]
    dyHx = dyu[1]
    dxHy = dxu[2]
    dxEz = dxu[3]
    dyEz = dyu[3]

    diff42_y!(dyEz, Ez, dy)
    diff42_x!(dxEz, Ez, dx)
    diff42_x!(dxHy, Hy, dy)
    diff42_y!(dyHx, Hx, dx)

    @. dtHx = -dyEz
    @. dtHy =  dxEz
    @. dtEz =  dxHy - dyHx

    sommerfeld_bcs(dtu[1], u[1], x, y)
    sommerfeld_bcs(dtu[2], u[2], x, y)
    sommerfeld_bcs(dtu[3], u[3], x, y)
end

function dirichlet_bcs(dtu, u, x, y)
    nx, ny = size(u)

    ############  j = 1
    j = 1
    for i = 1:nx
        dtu[i,j] = 0.0
    end

    ############  j = ny
    j = ny
    for i = 1:nx
        dtu[i,j] = 0.0
    end

    ############  i = 1
    i = 1
    for j = 1:ny
        dtu[i,j] = 0.0
    end

    ############  i = nx
    i = nx
    for j = 1:ny
        dtu[i,j] = 0.0
    end
 
end

function sommerfeld_bcs_deriv(dtu, u, dxu, dyu, x, y)

    nx, ny = size(u)

    u0::Float64 = 0.0
    u_falloff::Float64 = 2.0
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    idx2 = 1.0/(2.0*dx)
    idy2 = 1.0/(2.0*dy)

    diff22_x!(dxu, u, dx)
    diff22_y!(dyu, u, dy)

    ############  j = 1
    j = 1
    for i = 1:nx
        dtu[i,j] = - (x[i]*dxu[i,j] + y[j]*dyu[i,j] + u_falloff*(u[i,j] - u0))/sqrt(x[i]*x[i] + y[j]*y[j])
    end

    ############  j = ny
    j = ny
    for i = 1:nx
        dtu[i,j] = - (x[i]*dxu[i,j] + y[j]*dyu[i,j] + u_falloff*(u[i,j] - u0))/sqrt(x[i]*x[i] + y[j]*y[j])
    end

    ############  i = 1
    i = 1
    for j = 1:ny
        dtu[i,j] = - (x[i]*dxu[i,j] + y[j]*dyu[i,j] + u_falloff*(u[i,j] - u0))/sqrt(x[i]*x[i] + y[j]*y[j])
    end

    ############  i = nx
    i = nx
    for j = 1:ny
        dtu[i,j] = - (x[i]*dxu[i,j] + y[j]*dyu[i,j] + u_falloff*(u[i,j] - u0))/sqrt(x[i]*x[i] + y[j]*y[j])
    end
    
end

function sommerfeld_bcs(dtu, u, x, y)

    nx, ny = size(u)

    u0::Float64 = 0.0
    u_falloff::Float64 = 2.0
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    idx2 = 1.0/(2.0*dx)
    idy2 = 1.0/(2.0*dy)

    #########################  j=1 
    j = 1

    i = 1
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx2
    dyu = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    for i = 2:nx-1
        inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
        dxu = (u[i+1,j] - u[i-1,j])*idx2
        dyu = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy2
        dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))
    end

    i = nx
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = ( 3.0*u[nx,j] - 4.0*u[nx-1,j] + u[nx-2,j]) * idx2
    dyu = (-3.0*u[i,1]  + 4.0*u[i,2]    - u[i,3]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    #########################  j=ny
    j = ny

    i = 1
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx2
    dyu = ( 3.0*u[i,ny] - 4.0*u[i,ny-1] + u[i,ny-2]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    for i = 2:nx-1
        inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
        dxu = (u[i+1,j] - u[i-1,j])*idx2
        dyu = ( 3.0*u[i,ny] - 4.0*u[i,ny-1] + u[i,ny-2]) * idy2
        dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))
    end

    i = nx
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = ( 3.0*u[nx,j] - 4.0*u[nx-1,j] + u[nx-2,j]) * idx2
    dyu = ( 3.0*u[i,ny] - 4.0*u[i,ny-1] + u[i,ny-2]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    #########################  i=1 
    i = 1

    j = 1
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx2
    dyu = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    for j = 2:ny-1
        inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
        dxu = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx2
        dyu = (u[i,j+1] - u[i,j-1])*idy2
        dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))
    end

    j = ny
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = (-3.0*u[1,j] + 4.0*u[2,j] - u[3,j]) * idx2
    dyu = ( 3.0*u[i,ny] - 4.0*u[i,ny-1] + u[i,ny-2]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    #########################  i=nx
    i = nx

    j = 1
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = ( 3.0*u[nx,j] - 4.0*u[nx-1,j] + u[nx-2,j]) * idx2
    dyu = (-3.0*u[i,1] + 4.0*u[i,2] - u[i,3]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

    for j = 2:ny-1
        inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
        dxu = ( 3.0*u[nx,j] - 4.0*u[nx-1,j] + u[nx-2,j]) * idx2
        dyu = (u[i,j+1] - u[i,j-1])*idy2
        dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))
    end

    j = ny
    inv_r = 1.0/sqrt(x[i]*x[i] + y[j]*y[j])
    dxu = ( 3.0*u[nx,j] - 4.0*u[nx-1,j] + u[nx-2,j]) * idx2
    dyu = ( 3.0*u[i,ny] - 4.0*u[i,ny-1] + u[i,ny-2]) * idy2
    dtu[i,j] = - inv_r*(x[i]*dxu + y[j]*dyu + u_falloff*(u[i,j] - u0))

end
