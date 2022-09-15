module NBody

using DifferentialEquations

export Bodies, Params, init_data_binary, init_data_3bodies, intNbdyEqs

struct Bodies

    ndim :: Int64
    nbdy :: Int64
    mass :: Vector{Float64}
    u :: Vector{Float64}
    G :: Float64

    function Bodies(q0::Vector{Vector{Float64}}, p0, mass, G)
        nbdy = length(mass)
        ndim = length(q0[1])

        idx::Int64 = 0
        u = zeros(Float64,2*nbdy*ndim)
        for x in q0
            for y in x
                idx += 1
                u[idx] = y
            end
        end
        for x in p0
            for y in x
                idx += 1
                u[idx] = y
            end
        end
        new(ndim, nbdy, mass, u, G)
    end
end

struct Params
    tf   :: Float64
    abin :: Float64
    ebin :: Float64
    a3   :: Float64
    e3   :: Float64
    function Params(tf)
        new(tf, 2.0, 0.0, 4.0, 0.5)
    end
    function Params(tf, abin, ebin, a3, e3)
        new(tf, abin, ebin, a3, e3)
    end
end

function init_data_binary(ndim, m1, m2, a, e, G=1.0)
    M = m1 + m2
    r = a*(1.0 - e)
    xa1 = m2/M*r
    xb1 = -m1/M*r
    pax = 0.0
    pay = m1*m2/M*sqrt(G*M*(1.0+e)/(a*(1.0-e)))
    pbx = -pax
    pby = -pay
    period = 2.0*pi*sqrt(a^3/(G*M))

    if ndim == 2
        qdata = [Float64[xa1, 0.0], Float64[xb1, 0.0]]
        pdata = [Float64[pax, pay], Float64[pbx, pby]]
    elseif ndim ==3
        qdata = [Float64[xa1, 0.0, 0.0], Float64[xb1, 0.0, 0.0]]
        pdata = [Float64[pax, pay, 0.0], Float64[pbx, pby, 0.0]]
    else
        throw(AssertionError("Require: ndim=2 or ndim=3"))
    end

    return qdata, pdata, period
end

function init_data_3bodies(ndim, m1, m2, m3, a1, e1, a3, e3, G=1.0)

    qbin, pbin, taubin = init_data_binary(ndim, m1, m2, a1, e1, G)
    q1 = qbin[1]
    q2 = qbin[2]
    p1 = pbin[1]
    p2 = pbin[2]

    M = m1 + m2
    qcm = (m1*q1+m2*q2)/M
    pcm = p1 + p2

    qbin3, pbin3, taubin3 = init_data_binary(ndim, M, m3, a3, e3, G)
    qb = qbin3[1]
    q3 = qbin3[2]
    pb = qbin3[1]
    p3 = pbin3[2]

    q1 += qb
    q2 += qb
    p1 += pb
    p2 += pb

    M3 = m1 + m2 + m3
    qcm = (m1*q1 + m2*q2 + m3*q3)/M3
    vcm = (p1 + p2 + p3)/M3
    q1 -= qcm
    q2 -= qcm
    q3 -= qcm

    p1 -= m1*vcm
    p2 -= m2*vcm
    p3 -= m3*vcm
    return [q1, q2, q3],[p1, p2, p3], max(taubin,taubin3)
end

function norm(u)
    s = zero(eltype(u))
    for x in u
        s += x*x
    end
    return sqrt(s)
end

function nbdy_eqs!(dtu, u, b, t)

    # epsilon is the softening length
    epsq::Float64 = (1.0e-4)^2
    
    ndim = b.ndim
    m = b.mass
    nbdy = length(m)
    G = b.G

    q = zeros(Float64,ndim,nbdy)
    p = zeros(Float64,ndim,nbdy)
    dtq = zeros(Float64,ndim,nbdy)
    dtp = zeros(Float64,ndim,nbdy)

    r = zeros(Float64,ndim)
    idx::Int64 = 0
    for i = 1:nbdy
        for j = 1:ndim
            idx += 1
            q[j,i] = u[idx]
        end
    end
    for i = 1:nbdy
        for j = 1:ndim
            idx += 1
            p[j,i] = u[idx]
        end
    end

    for i = 1:nbdy
        dtq[:,i] = p[:,i]/m[i]
        for j = i+1:nbdy
            r[:] = q[:,j] - q[:,i]
            rsq = sum(r.*r) + epsq
            dtp[:,i] += G*m[i]*m[j]*r[:]/rsq^1.5
            dtp[:,j] -= dtp[:,i]
        end
    end

    #=
    println("position 1: ",q[1,1],"  ",q[2,1])
    println("position 2: ",q[1,2],"  ",q[2,2])
    println("momentum 1: ",p[1,1],"  ",p[2,1])
    println("r         : ",r[1],"  ",r[2],"  ",norm(r))
    =#
    
    
    idx = 0
    for i = 1:nbdy
        for j = 1:ndim
            idx += 1
            dtu[idx] = dtq[j,i]
        end
    end
    for i = 1:nbdy
        for j = 1:ndim
            idx += 1
            dtu[idx] = dtp[j,i]
        end
    end
end

function intNbdyEqs(tspan, dt, b::Bodies)
   
    u = b.u

    prob = ODEProblem(nbdy_eqs!, u, tspan, b)
    sol = solve(prob, Vern9(), reltol=1.0e-12, abstol=1.0e-12, saveat=dt);

    return sol
end


end #module
