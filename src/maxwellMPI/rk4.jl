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

    gh = fields.gh
    dt = gh.dt
    xi = gh.lcoords
    dxi = gh.dx0

    un = fields.u
    utmp = fields.u2
    k1 = fields.wrk[1]
    k2 = fields.wrk[2]
    k3 = fields.wrk[3]
    k4 = fields.wrk[4]
    dxu = fields.dxu
    dyu = fields.dyu
    comm = gh.comm

    dthalf = 0.5*dt
    thalf = tx + dthalf
    tnp1 = tx + dt
    dvars = fields.dvars

    func(k1, un, dxu, dyu, xi, dxi, dvars, tx)
    @. rk4_helper(utmp, un, k1, dthalf)
    waveguide_bcs(utmp)
    Maxwell2D.grid_sync!(utmp, gh, comm)

    func(k2, utmp, dxu, dyu, xi, dxi, dvars, thalf)
    @. rk4_helper(utmp, un, k2, dthalf)
    waveguide_bcs(utmp)
    Maxwell2D.grid_sync!(utmp, gh, comm)

    func(k3, utmp, dxu, dyu, xi, dxi, dvars, thalf)
    @. rk4_helper(utmp, un, k3, dt)
    waveguide_bcs(utmp)
    Maxwell2D.grid_sync!(utmp, gh, comm)

    func(k4, utmp, dxu, dyu, xi, dxi, dvars, tnp1)
    @. rk4_helper2(un, k1, k2, k3, k4, dt)
    waveguide_bcs(un)
    Maxwell2D.grid_sync!(un, gh, comm)

    t[1] = tnp1

end

