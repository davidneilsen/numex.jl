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
    dz = fields.grid.dz
    dt = fields.grid.dt
    x  = fields.grid.x
    y  = fields.grid.y
    z  = fields.grid.z
    yn = fields.u
    ynp1 = fields.u2
    f = fields.wrk
    dxy = fields.dxu
    dyy = fields.dyu
    dzy = fields.dzu
    
    func(f, yn, dxy, dyy, dzy, x, y, z, dx, dy, dz, tx)
    @. kodiss!(f, yn, dx, dy, dz)
    #@printf("h. f1=%g, f2=%g\n",l2norm(f1[1]),l2norm(f1[2]))
    @. rk2_helper1(ynp1, yn, f, dt)
    
    thalf = tx + 0.5*dt
    t[1] += dt
    func(f, ynp1, dxy, dyy, dzy, x, y, z, dx, dy, dz, thalf)
    @. kodiss!(f, yn, dx, dy, dz)
    @. rk2_helper2(yn, f, dt)
end

