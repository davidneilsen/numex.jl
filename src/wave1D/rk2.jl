#=------------------------------------------------------------------
 =
 =  The classic 2nd Order Runge-Kutta integrator, or the Midpoint method
 =
 =  The PDES are solved using the Method of Lines.  The variables defined
 =  in fields.u are integrated in time.  fields.u2 is used as work space.
 =  The solution at the advanced time is written back into fields.u,
 =  so fields.u is changed by this routine.
 =
 =------------------------------------------------------------------=#

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
    dt = fields.grid.dt
    x  = fields.grid.x
    yn = fields.u
    ynp1 = fields.u2
    f = fields.wrk
    dxy = fields.dxu
    
    func(f, yn, dxy, x, dx, tx)
    @. kodiss!(f, yn, dx)
    @. rk2_helper1(ynp1, yn, f, dt)
    
    thalf = tx + 0.5*dt
    t[1] += dt
    func(f, ynp1, dxy, x, dx, thalf)
    @. kodiss!(f, yn, dx)
    @. rk2_helper2(yn, f, dt)
end

