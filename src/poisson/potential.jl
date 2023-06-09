#-----------------------------------------------------------------------
#  Parameter structure.
#
#  shp = [nx, ny]  number of points in each dimension
#  bbox = [xmin, xmax, ymin, ymax]  grid domain limits
#  dielectrics = [d1, d2, ... ]  array of dielectric objects.  Currently
#                                these can be rectangles or circles.
#                                Each data structure is an 8-element real array.
#  Vb = [Vtop, Vright, Vbottom, Vleft]  Dirichlet boundary conditions for
#                                       a rectangular domain.
#  V0 = initial guess for potential
#  perm1 = default permittivity (vacuum value)  Set to 1.0
#  MASK_LIVE = mask value for points that should be updated
#  MASK_BOUND = mask value for points with fixed potential.  Not updated.
#
#-----------------------------------------------------------------------
struct Pars
    shp::Vector{Int64}
    bbox::Vector{Float64}
    
    dielectrics::Vector{Vector{Float64}}
    
    Vb::Vector{Float64}
    
    V0::Float64
    perm1::Float64
    
    MASK_LIVE::Int64
    MASK_BOUND::Int64
    
    function Pars(shp, bbox, d1, Vb)
        V0 = 1.0
        p1 = 1.0        
        new(shp, bbox, d1, Vb, V0, p1, 1, 0)
    end
end

#-----------------------------------------------------------------------
#  Sets the permittivity for a dielectric rectangular region.
#  Not called by the user.
#  pars = dielectric data structure, not the general par structure
#-----------------------------------------------------------------------
function dielectricSquare!(perm, x, y, pars)
    nx, ny = size(perm)
    
    x1 = pars[2]
    x2 = pars[3]
    y1 = pars[4]
    y2 = pars[5]
    p1 = pars[6]
    
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    
    i1 = trunc(Int,(x1 - x[1])/dx) + 1
    i2 = trunc(Int,(x2 - x[1])/dx) + 1
    j1 = trunc(Int,(y1 - y[1])/dy) + 1
    j2 = trunc(Int,(y2 - y[1])/dy) + 1    
    
    for j = j1:j2
        for i = i1:i2
            perm[i,j] = p1
        end
    end
end


#-----------------------------------------------------------------------
#  Sets the permittivity for a dielectric cicular region.
#  Not called by the user.
#  pars = dielectric data structure, not the general par structure
#-----------------------------------------------------------------------
function dielectricCyl!(perm, x, y, pars)
    nx, ny = size(perm)
    
    cx0 = pars[2]
    cy0 = pars[3]
    crad = pars[4]
    p1 = pars[5]
    
    crsq = crad^2
    
    for j = 1:ny
        for i = 1:nx
            r2 = (x[i] - cx0)^2 + (y[j] - cy0)^2
            if r2 < crsq
                perm[i,j] = p1
            end
        end
    end
end

#-----------------------------------------------------------------------
#  Initialize the potential, permittivity, and the mask.  This is the
#  primary routine for initialization.  Called by the solver routine.
#  It sets the Dirichlet BCs on the outer boundaries, and calls routines
#  to set up the dielectrics.
#-----------------------------------------------------------------------
function initPotentialDielectric!(V, perm, mask, x, y, pars)
    nx, ny = size(V)
    
    V0 = pars.V0
    p1 = pars.perm1
    MASK_BOUND = pars.MASK_BOUND
    MASK_LIVE = pars.MASK_LIVE

    Vtop = pars.Vb[1]
    Vbottom = pars.Vb[3]
    Vleft = pars.Vb[4]
    Vright = pars.Vb[2]
    
    dx = x[2] - x[1]
    dy = y[2] - y[1] 
    
    # default values everywhere
    for j = 1:ny
        for i = 1:nx
            V[i,j] = V0
            mask[i,j] = MASK_LIVE
            perm[i,j] = p1
        end
    end

    # bottom boundary
    for i = 1:nx
        V[i,1] = Vbottom
        mask[i,1] = MASK_BOUND
    end
    
    # top boundary
    for i = 1:nx
        V[i,ny] = Vtop
        mask[i,ny] = MASK_BOUND
    end
    
    # left boundary
    for j = 1:ny
        V[1,j] = Vleft
        mask[1,j] = MASK_BOUND
    end
    
    # right boundary
    for j = 1:ny
        V[nx,j] = Vright
        mask[nx,j] = MASK_BOUND
    end
    
    # dielectric
    d1 = pars.dielectrics
    for m = 1:length(d1)
        dielec = d1[m]
        if dielec[1] > 0.9 && dielec[1] < 1.1
            dielectricSquare!(perm, x, y, dielec)
        elseif dielec[1] > 1.9 && dielec[1] < 2.1
            dielectricCyl!(perm, x, y, dielec)
        else
            println("Unknown dielectric type")
        end
    end      
end

#-----------------------------------------------------------------------
#  Initialize the permittivity and mask.  This is the
#  primary routine for initialization.  Called by the solver routine.
#  It sets the Dirichlet BCs on the outer boundaries, and calls routines
#  to set up the dielectrics.
#-----------------------------------------------------------------------
function initDielectric!(perm, mask, x, y, pars)
    nx, ny = size(perm)
    
    p1 = pars.perm1
    MASK_BOUND = pars.MASK_BOUND
    MASK_LIVE = pars.MASK_LIVE

    dx = x[2] - x[1]
    dy = y[2] - y[1] 
    
    # default values everywhere
    for j = 1:ny
        for i = 1:nx
            mask[i,j] = MASK_LIVE
            perm[i,j] = p1
        end
    end

    # dielectric
    d1 = pars.dielectrics
    for m = 1:length(d1)
        dielec = d1[m]
        if dielec[1] > 0.9 && dielec[1] < 1.1
            dielectricSquare!(perm, x, y, dielec)
        elseif dielec[1] > 1.9 && dielec[1] < 2.1
            dielectricCyl!(perm, x, y, dielec)
        else
            println("Unknown dielectric type")
        end
    end      
end


#-----------------------------------------------------------------------
#
#
#-----------------------------------------------------------------------
function calE!(Ex, Ey, V, dx, dy; rescale=false, scale=1.0)
    deriv_x21(Ex, V, dx)
    deriv_y21(Ey, V, dy)
    @. Ex = -Ex
    @. Ey = -Ey

    nx, ny = size(V)
    if rescale
        for j = 1:ny
            for i = 1:nx
                nrm = sqrt(Ex[i,j]^2 + Ey[i,j]^2) + 1.0e-10
                temp = scale/nrm
                Ex[i,j] *= temp
                Ey[i,j] *= temp
            end
        end
    end
end

#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
function setSqr(center, a, b, perm)
    dp = zeros(8)
    dp[1] = 1.0
    dp[2] = center[1] - 0.5*a
    dp[3] = center[1] + 0.5*a
    dp[4] = center[2] - 0.5*b
    dp[5] = center[2] + 0.5*b
    dp[6] = perm
    return dp
end

#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
function setCyl(center, rad, perm)
    dp = zeros(8)
    dp[1] = 2.0
    dp[2] = center[1]
    dp[3] = center[2]
    dp[4] = rad
    dp[5] = perm
    return dp
end

#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
function setCharge(center, charge)
    dp = zeros(4)
    dp[1] = 1.0  # point charge
    dp[2] = center[1]
    dp[3] = center[2]
    dp[4] = charge
    return dp
end

#-----------------------------------------------------------------------
#  This routine solves the Laplace equation with Successive Over Relaxation.
#  It also performs the update with a Red-Black ordering of grid points.
#  It is much faster than Gauss-Seidel iteration.
#-----------------------------------------------------------------------
function SOR!(V,omega,p,mask,tol,maxSteps,h,pars)

    nx, ny = size(V)
    step = 0
    
    err = 2.0*tol
    MASK_LIVE = pars.MASK_LIVE
    
    while (step < maxSteps && err > tol)
        #print("begin loop, tol=",tol)
        errsum = 0.0
        for phase = 0:1
            for j = 2:ny-1
                ibgn = 1 + mod(j+phase,2)
                for i = ibgn:2:nx-1
                    Res = 0.0
                    if mask[i,j] == MASK_LIVE
                        a0 = 1.0/(p[i,j] + p[i-1,j] + p[i,j-1] + p[i-1,j-1])
                        a1 = 0.5*(p[i,j] + p[i,j-1])
                        a2 = 0.5*(p[i-1,j] + p[i,j])
                        a3 = 0.5*(p[i-1,j-1] + p[i-1,j])
                        a4 = 0.5*(p[i,j-1] + p[i-1,j-1])
                    
                        Res = a0*(a1*V[i+1,j] + a2*V[i,j+1] + a3*V[i-1,j] + a4*V[i,j-1]) - V[i,j]
                    end
                    errsum += (omega*Res)^2
                    V[i,j] += omega*Res
                end
            end
        end
        err = errsum*h^2
        step += 1
    end
 
   return step, err
end

#-----------------------------------------------------------------------
#
#
#-----------------------------------------------------------------------
function defA(eps, Vbound)
    nx, ny = size(eps)
    
    Vtop = Vbound[1]
    Vright = Vbound[2]
    Vbottom = Vbound[3]
    Vleft = Vbound[4]
    
    nd = nx*ny
    A = zeros(Float64,nd,nd)
    f = zeros(Float64,nd)
    
    # add a point charge as a test.  Need to generalize
    #nmid = trunc(Int64,nx/2)
    #f[nmid + (nmid-1)*nx] = -109.0
    
    for j = 1:ny
        for i = 1:nx
            ii = i + (j-1)*nx
            ip = (i+1) + (j-1)*nx
            im = (i-1) + (j-1)*nx
            jp = i + j*nx
            jm = i + (j-2)*nx
            
            eim = 1.0
            ejp = 1.0
            ejm = 1.0
            ejp = 1.0
            
            if i == 1 
                # left boundary
                eim = eps[i,j]
                eip = 0.5*(eps[i,j] + eps[i+1,j])
                f[ii] -= eim*Vleft
                A[ii,ip] = eip
            elseif i == nx
                # right boundary
                eip = eps[i,j]
                eim = 0.5*(eps[i,j] + eps[i-1,j])
                f[ii] -= eip*Vright
                A[ii,im] = eim
            else
                eip = 0.5*(eps[i,j] + eps[i+1,j])
                eim = 0.5*(eps[i,j] + eps[i-1,j])
                A[ii,ip] = eip
                A[ii,im] = eim
            end
            
            if j == 1
                # bottom boundary
                ejm = eps[i,j]
                ejp = 0.5*(eps[i,j] + eps[i,j+1])
                f[ii] -= ejm*Vbottom
                A[ii,jp] = ejp
            elseif j == ny
                # top boundary
                ejp = eps[i,j]
                ejm = 0.5*(eps[i,j] + eps[i,j-1])
                f[ii] -= ejp*Vtop
                A[ii,jm] = ejm
            else
                # center
                ejm = 0.5*(eps[i,j] + eps[i,j-1])
                ejp = 0.5*(eps[i,j] + eps[i,j+1])
                A[ii,jm] = ejm
                A[ii,jp] = ejp
            end
            a0 = -(eip + eim + ejp + ejm)
            A[ii,ii] = a0
        end
    end
    return A, f
end
