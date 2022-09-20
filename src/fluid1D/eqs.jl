function init_p!(P, rho, pars)
    k = pars.kappa
    gamma = pars.gamma
    @. P = k * rho^gamma
end

function p2c!(u, v, x, pars)
    nx = length(u)
    D = u[1]
    S = u[2]
    tau = u[3]
    rho = v[1]
    vx = v[2]
    P = v[3]
    gamma = pars.gamma
    for i = 1:nx
        e = rho[i] + P[i]/(gamma - 1.0)
        Wsq = 1.0/(1.0 - vx[i]^2)
        D[i] = rho[i]*sqrt(Wsq)
        S[i] = (e + P[i]) * Wsq * vx[i]
        tau[i] = (e + P[i]) * Wsq - P[i] - D[i]
    end
end

function set_initial!(fields)
    v = fields.v
    nx = length(v[1])
    m::Int64 = floor(nx/2)

    rho  = v[1]
    vx = v[2]
    P = v[3]
    for i = 1:m
        rho[i] = id.rho_left
        vx[i] = id.vx_left
        P[i] = id.P_left
    end
    for i = m+1:nx
        rho[i] = id.rho_right
        vx[i] = id.vx_right
        P[i] = id.P_right
    end

    u = fields.u
    x = fields.grid.x
    pars = fields.pars
    p2c!(u, v, x, pars)
end
