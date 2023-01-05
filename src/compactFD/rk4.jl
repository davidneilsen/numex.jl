function rk4_helper(y2, y1, f, dt)
    @. y2 = y1 + dt*f
end

function rk4_helper2(un, f1, f2, f3, f4, dt)
    t1 = dt / 6.0
    @. un += t1*(f1 + 2.0*f2 + 2.0*f3 + f4)
end

function rk4_step!(func::Function, fields, t)
    # func calculates the RHS.  The boundary conditions are returned in func.
    neqs = length(fields.u)
    tx = t[1]

    dx = fields.grid.dx
    dt = fields.grid.dt

    x  = fields.grid.x
    un = fields.u
    utmp = fields.ut
    f1 = fields.k1
    f2 = fields.k2
    f3 = fields.k3
    f4 = fields.k4
    dxu = fields.dxu
    A = fields.A
    B = fields.B
    F1 = fields.F1
    F2 = fields.F2
    order = fields.order
    method = fields.method
    
    dthalf = 0.5*dt
    thalf = tx + dthalf
    tnp1 = tx + dt
    t[1] = tnp1

    func(f1, un, dxu, x, dx, tx, A, B, order, method)

    @. rk4_helper(utmp, un, f1, dthalf)
    func(f2, utmp, dxu, x, dx, thalf, A, B, order, method)

    @. rk4_helper(utmp, un, f2, dthalf)
    func(f3, utmp, dxu, x, dx, thalf, A, B, order, method)

    @. rk4_helper(utmp, un, f3, dt)
    func(f4, utmp, dxu, x, dx, tnp1, A, B, order, method)

    @. rk4_helper2(un, f1, f2, f3, f4, dt)

end

