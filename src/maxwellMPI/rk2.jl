function rk2_helper1(ynp1, yn, f, dt)
    @. ynp1 = yn + 0.5*dt*f
end
function rk2_helper2(yn, f, dt)
    @. yn += dt[1]*f
end

function rk2_step!(func::Function, fields, t)
    # func calculates the RHS.  The boundary conditions are returned in func.
    neqs = length(fields.u)
    tx = t[1]

    dx = fields.grid.dx
    dy = fields.grid.dy
    dt = fields.grid.dt
    x  = fields.grid.x
    y  = fields.grid.y
    yn = fields.u
    ynp1 = fields.u2
    f = fields.wrk
    dxy = fields.dxu
    dyy = fields.dyu
    
    func(f, yn, dxy, dyy, x, y, dx, dy, tx)
    @. kodiss!(f, yn, dx, dy)
    #@printf("h. f1=%g, f2=%g\n",l2norm(f1[1]),l2norm(f1[2]))
    @. rk2_helper1(ynp1, yn, f, dt)
    Maxwell2D.grid_sync!(fields, comm)
    
    thalf = tx + 0.5*dt
    t[1] += dt
    func(f, ynp1, dxy, dyy, x, y, dx, dy, thalf)
    @. kodiss!(f, yn, dx, dy)
    @. rk2_helper2(yn, f, dt)
    Maxwell2D.grid_sync!(fields, comm)
end

function rk2_step_dummy!(func::Function, fields, t)
    dt = fields.gh.dt
    dx = fields.gh.dx0[1]
    dy = fields.gh.dx0[2]
    un = fields.u
    f = fields.wrk
    comm = fields.gh.comm

    t[1] += dt
    dxu = fields.dxu
    dyu = fields.dyu

    nx, ny = size(un[1])
    for m = 1:3
        for j = 1:ny
            for i = 1:nx
                f[m][i,j] = un[m][i,j]
            end
        end
    end
    for i = 1:3
        diff22_x!(un[i], f[i], dy)
        @printf("rk2_step_dumm! Pre-sync >>> i=%d, |un[i]|=%g\n",i,l2norm(un[i]))
    end
    
    Maxwell2D.grid_sync!(fields, comm)

    for i = 1:3
        @printf("rk2_step_dumm! Post-sync >>> i=%d, |un[i]|=%g\n",i,l2norm(un[i]))
    end
    
end
