using LinearAlgebra

function init_data!(fields, t, pars)

    idtype = pars["idtype"]
    if idtype == 0
        gauss_init_data!(fields)
    elseif idtype == 1
        a = fields.gh.cbox0[1,2]
        b = fields.gh.cbox0[2,2]
        x = fields.gh.lcoords[1]
        y = fields.gh.lcoords[2]
        waveguide_init_data!(fields.u, x, y, t, a, b, pars)

        proc = fields.proc
        rank = fields.gh.rank
        @. proc = rank
    else
        Maxwell2D.oops("unknown idtype")
    end

end

function cal_error!(fields, t, pars)
    idtype = pars["idtype"]

    if idtype == 0
        for i = 1:m
            @. fields.u2[m] = 0.0;
        end
    elseif idtype == 1
        a = fields.gh.cbox0[1,2]
        b = fields.gh.cbox0[2,2]
        x = fields.gh.lcoords[1]
        y = fields.gh.lcoords[2]
        waveguide_init_data!(fields.dxu, x, y, t, a, b, pars)
        for m = 1:3
            @. fields.u2[m] = fields.u[m] - fields.dxu[m]
        end
    else
        Maxwell2D.oops("unknown idtype")
    end

end

function Xwaveguide_init_data!(fields)

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]
    x = fields.gh.lcoords[1]
    y = fields.gh.lcoords[2]
    shp = size(Ex)
    time = 0.0

    for j=1:shp[2], i=1:shp[1]
        Ex[i,j] = -cos(sqrt(2.0)*pi*time)*cos(pi*x[i])*sin(pi*y[j])
        Ey[i,j] =  cos(sqrt(2.0)*pi*time)*sin(pi*x[i])*cos(pi*y[j])
        Bz[i,j] =  -sqrt(2.0)*sin(sqrt(2.0)*pi*time)*cos(pi*x[i])*cos(pi*y[j])
    end

    proc = fields.proc
    rank = fields.gh.rank
    @. proc = rank
end

function waveguide_init_data!(u, x, y, t, a, b, pars)

    Ex = u[1]
    Ey = u[2]
    Bz = u[3]
    shp = size(Ex)

    B0 = pars["idwg_bzamp"]
    m = pars["idwg_m_mode"]
    n = pars["idwg_n_mode"]
    omn = pi*sqrt((m/a)^2 + (n/b)^2)

    for j=1:shp[2], i=1:shp[1]
        Ex[i,j] = -B0*n*pi*cos((m*pi*x[i])/a)*sin(omn*t)*sin((n*pi*y[j])/b)/(b*omn)
        Ey[i,j] =  B0*m*pi*cos((n*pi*y[j])/b)*sin(omn*t)*sin((m*pi*x[i])/a)/(a*omn)
        Bz[i,j] =  B0*cos((m*pi*x[i])/a)*cos((n*pi*y[j])/b)*cos(omn*t)
    end
end


function gauss_init_data!(fields)
    nx = fields.gh.lshp[1]
    ny = fields.gh.lshp[2]

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]
    x = fields.gh.lcoords[1]
    y = fields.gh.lcoords[2]
    proc = fields.proc
    rank = fields.gh.rank

    amp1 = 1.0
    lambda1 = 0.4

    for j=1:ny, i=1:nx
        r = sqrt(x[i]*x[i] + y[j]*y[j])
        Bz[i,j] = - 8.0*amp1*lambda1*lambda1*exp(-r*r*lambda1)
    end
    for j=1:ny, i=1:nx
        proc[i,j] = rank
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

function maxwell_TE!(dtu, u, dxu, dyu, xi, dxi, dvars, time)
    x = xi[1]
    y = xi[2]
    dx = dxi[1]
    dy = dxi[2]

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
    shp = size(Ex)
    nx = shp[1]
    ny = shp[2]

    dtype = dvars.dtype
    Dx1d = dvars.Dx1d
    Dy1d = dvars.Dy1d
    Dx2d = dvars.Dx2d
    Dy2d = dvars.Dy2d
    u1x = dvars.u1x
    u1y = dvars.u1y
    du1x = dvars.du1x
    du1y = dvars.du1y

    if dvars.dissipation == -2
        # Apply compact filter to vars before taking derivatives
        for m = 1:3
            cfilter!(u[m], dvars)
        end
    end

    if dtype == 0
        diff22_x!(dxHz, Hz, dx)
        diff22_y!(dyHz, Hz, dy)
        diff22_x!(dxEy, Ey, dx)
        diff22_y!(dyEx, Ex, dy)
    elseif dtype == 1
        diff42_x!(dxHz, Hz, dx)
        diff42_y!(dyHz, Hz, dy)
        diff42_x!(dxEy, Ey, dx)
        diff42_y!(dyEx, Ex, dy)
    elseif dtype == 2
        diff666_x!(dxHz, Hz, dx)
        diff666_y!(dyHz, Hz, dy)
        diff666_x!(dxEy, Ey, dx)
        diff666_y!(dyEx, Ex, dy)
    elseif dtype == 3
        Hz1d = vec(Hz)
        Ey1d = vec(Ey)
        Ex1d = vec(Ex)
        dxHz1d = Dx2d*Hz1d
        dyHz1d = Dy2d*Hz1d
        dxEy1d = Dx2d*Ey1d
        dyEx1d = Dy2d*Ex1d
        dxEy = reshape(dxEy1d, nx, ny)
        dyEx = reshape(dyEx1d, nx, ny)
        dxHz = reshape(dxHz1d, nx, ny)
        dyHz = reshape(dyHz1d, nx, ny)
    elseif dtype == 4
        cdiff_x!(dxEy, Ey, Dx1d, u1x, du1x)
        cdiff_y!(dyEx, Ex, Dy1d, u1y, du1y)
        cdiff_x!(dxHz, Hz, Dx1d, u1x, du1x)
        cdiff_y!(dyHz, Hz, Dy1d, u1y, du1y)
    else
        
    end

    @. dtEx = dyHz
    @. dtEy = -dxHz
    @. dtHz = dyEx - dxEy

#
    # last using these
    #sommerfeld_bcs(dtu[1], u[1], x, y)
    #sommerfeld_bcs(dtu[2], u[2], x, y)
    #sommerfeld_bcs(dtu[3], u[3], x, y)
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

# 
function waveguide_bcs(u)
    Ex = u[1]
    Ey = u[2]
    Bz = u[3]
    shp = size(Ex)
    nx = shp[1]
    ny = shp[2]

    j = 1
    for i = 1:nx
        Ex[i,j] = 0.0
        Ey[i,j] = (48.0*Ey[i,j+1] - 36.0*Ey[i,j+2] + 16.0*Ey[i,j+3] - 3.0*Ey[i,j+4])/25.0
        Bz[i,j] = (48.0*Bz[i,j+1] - 36.0*Bz[i,j+2] + 16.0*Bz[i,j+3] - 3.0*Bz[i,j+4])/25.0
    end
 
    j = ny
    for i = 1:nx
        Ex[i,j] = 0.0
        Ey[i,j] = (48.0*Ey[i,j-1] - 36.0*Ey[i,j-2] + 16.0*Ey[i,j-3] - 3.0*Ey[i,j-4])/25.0
        Bz[i,j] = (48.0*Bz[i,j-1] - 36.0*Bz[i,j-2] + 16.0*Bz[i,j-3] - 3.0*Bz[i,j-4])/25.0
    end

    i = 1
    for j = 1:ny
        Ey[i,j] = 0.0
        Ex[i,j] = (48.0*Ex[i+1,j] - 36.0*Ex[i+2,j] + 16.0*Ex[i+3,j] - 3.0*Ex[i+4,j])/25.0
        Bz[i,j] = (48.0*Bz[i+1,j] - 36.0*Bz[i+2,j] + 16.0*Bz[i+3,j] - 3.0*Bz[i+4,j])/25.0
    end

    i = nx
    for j = 1:ny
        Ey[i,j] = 0.0
        Ex[i,j] = (48.0*Ex[i-1,j] - 36.0*Ex[i-2,j] + 16.0*Ex[i-3,j] - 3.0*Ex[i-4,j])/25.0
        Bz[i,j] = (48.0*Bz[i-1,j] - 36.0*Bz[i-2,j] + 16.0*Bz[i-3,j] - 3.0*Bz[i-4,j])/25.0
    end
end
