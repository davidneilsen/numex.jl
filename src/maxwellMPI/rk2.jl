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

    gh = fields.gh
    dx = gh.dx0[1]
    dy = gh.dx0[2]
    dt = gh.dt
    x  = gh.lcoords[1]
    y  = gh.lcoords[2]
    un = fields.u
    unp1 = fields.u2
    f = fields.wrk[1]
    dxu = fields.dxu
    dyu = fields.dyu
    comm = gh.comm
    dvars = fields.dvars

    xi = Vector{Vector{Float64}}(undef,2) 
    xi[1] = gh.lcoords[1]
    xi[2] = gh.lcoords[2]
    dxi = Vector{Float64}(undef,2)
    dxi[1] = gh.dx0[1]
    dxi[2] = gh.dx0[2]

    func(f, un, dxu, dyu, xi, dxi, dvars, tx)
    @. kodiss!(f, un, dx, dy)
    #@printf("h. f1=%g, f2=%g\n",l2norm(f1[1]),l2norm(f1[2]))
    @. rk2_helper1(unp1, un, f, dt)
    waveguide_bcs(unp1)
    Maxwell2D.grid_sync!(unp1, gh, comm)
    
    thalf = tx + 0.5*dt
    t[1] += dt
    func(f, unp1, dxu, dyu, xi, dxi, dvars, thalf)
    @. kodiss!(f, un, dx, dy)
    @. rk2_helper2(un, f, dt)
    waveguide_bcs(un)
    Maxwell2D.grid_sync!(un, gh, comm)
end

function rk2_step_dummy!(func::Function, fields, t)
    gh = fields.gh
    dt = gh.dt
    dx = gh.dx0[1]
    dy = gh.dx0[2]
    un = fields.u
    f = fields.wrk
    comm = gh.comm

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
    
    Maxwell2D.grid_sync!(un, gh, comm)

    for i = 1:3
        @printf("rk2_step_dumm! Post-sync >>> i=%d, |un[i]|=%g\n",i,l2norm(un[i]))
    end
    
end
