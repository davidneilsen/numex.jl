function rk4_helper(y2, y1, f, dt)
    @. y2 = y1 + dt*f
end

function rk4_helper2(un, f1, f2, f3, f4, dt)
    t1 = dt / 6.0
    @. un += t1*(f1 + 2.0*f2 + 2.0*f3 + f4)
end

function rk4_step!(func::Function, u, pars, t, dt)
    # func calculates the RHS.  The boundary conditions are returned in func.
    tx = t[1]

    utmp = similar(u)

    f1 = similar(u)
    f2 = similar(u)
    f3 = similar(u)
    f4 = similar(u)
    
    dthalf = 0.5*dt
    thalf = tx + dthalf
    tnp1 = tx + dt
    t[1] = tnp1

    func(f1, u, pars, tx)
    
    rk4_helper(utmp, u, f1, dthalf)
    func(f2, utmp, pars, thalf)
    
    rk4_helper(utmp, u, f2, dthalf)
    func(f3, utmp, pars, thalf)
    
    rk4_helper(utmp, u, f3, dt) 
    func(f4, utmp, pars, tnp1)
    
    rk4_helper2(u, f1, f2, f3, f4, dt)

end

#-------------------------------------RK2------------------------------------

function rk2_helper1(ynp1, yn, f, dt)
    @. ynp1 = yn + 0.5*dt*f
end
function rk2_helper2(yn, f, dt)
    @. yn += dt[1]*f
end

function rk2_step!(func::Function, yn, pars, t, dt)
    tx = t[1]
    
    ynp1 = similar(yn)
    f = similar(yn)
    
    func(f, yn, pars, tx)
    #@printf("h. f1=%g, f2=%g\n",l2norm(f1[1]),l2norm(f1[2]))
    rk2_helper1(ynp1, yn, f, dt)
    
    thalf = tx + 0.5*dt
    t[1] += dt
    func(f, ynp1, pars, thalf)
    rk2_helper2(yn, f, dt)
end 
