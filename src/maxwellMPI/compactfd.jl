#=
#===============================================
#
#
#
================================================#
struct Derivs

    Dx2d = Array

    function Derivs(dtype, dx0, shp)

        nx = shp[1]
        ny = shp[2]
        dx = dx0[1]
        dy = dx0[2]

        if dtype == 1
            Px, Qx = derivKim4(nx, dx)
            Py, Qy = derivKim4(ny, dy)
        elseif dtype == 2
        else
            @printf("unknown dtype = %d\n");
        end

        Dx1d = Px \ Qx
        Dy1d = Py \ Qy
        Dx2d = BandedBlockBandedMatrix(Kron(Eye(ny), Dx1d))
        Dy2d = BandedBlockBandedMatrix(Kron(Dy1d, Eye(nx)));
        new(Dx2d, Dy2d)
    end

end
=#

#===============================================
#
#
#
================================================#
function dT4(n, h)
    dl = fill(0.25,n-1)
    du = fill(0.25,n-1)
    d = ones(n)
    dl[n-1] = 3.0
    du[1] = 3.0
    A = Tridiagonal(dl, d, du)

    B = zeros(n,n)
    hm1 = 1.0/h
    t1 = 0.75*hm1
    for i = 2:n-1
        B[i-1,i] = t1
        B[i+1,i] = -t1
    end

    B[1,1] = -17.0/6.0*hm1
    B[2,1] = -t1
    B[1,2] = 1.5*hm1
    B[1,3] = 1.5*hm1
    B[1,4] = -1.0/6.0*hm1

    B[n-1,n] = t1
    B[n,n-3] = 1.0/6.0*hm1
    B[n,n-2] = -1.5*hm1
    B[n,n-1] = -1.5*hm1
    B[n,n] = 17.0/6.0*hm1

    return A, B
end

#===============================================
#
#
#
================================================#
function dT6(n, h)
    third = 1.0/3.0
    dl = fill(third, n-1)
    du = fill(third, n-1)
    d = ones(n)
    dl[1] = 2.0/11.0
    dl[n-2] = 2.0/11.0
    dl[n-1] = 5.0
    du[1] = 5.0
    du[2] = 2.0/11.0
    du[n-1] = 2.0/11.0
    A = Tridiagonal(dl, d, du)

    hm1 = 1.0/h
    B = zeros(n,n)
    t1 = 1.0/36.0*hm1
    t2 = 7.0/9.0*hm1
    for j = 3:n-2
        B[j-2,j] =  t1
        B[j-1,j] =  t2
        B[j+1,j] = -t2
        B[j+2,j] = -t1
    end

    B[1,1] = -197.0/60.0*hm1
    B[2,1] = -20.0/33.0*hm1
    B[3,1] = -t1

    B[1,2] = -5.0/12.0*hm1
    B[2,2] = -35.0/132.0*hm1
    B[3,2] = -t2
    B[4,2] = -t1

    B[1,3] = 5.0*hm1
    B[2,3] = 34.0/33.0*hm1

    B[1,4] = -5.0/3.0*hm1
    B[2,4] = -7.0/33.0*hm1

    B[1,5] = 5.0/12.0*hm1
    B[2,5] = 2.0/33.0*hm1

    B[1,6] = -1.0/20.0*hm1
    B[2,6] = -1.0/132.0*hm1

    B[n-1,n-5] = 1.0/132.0*hm1
    B[n,n-5] = 1.0/20.0*hm1

    B[n-1,n-4] = -2.0/33.0*hm1
    B[n,n-4] = -5.0/12.0*hm1

    B[n-1,n-3] = 7.0/33.0*hm1
    B[n,n-3] = 5.0/3.0*hm1

    B[n-1,n-2] = -34.0/33.0*hm1
    B[n,n-2] = -5.0*hm1

    B[n-3,n-1] = t1
    B[n-2,n-1] = t2
    B[n-1,n-1] = 35.0/132.0*hm1
    B[n,n-1] = 5.0/12.0*hm1

    B[n-2,n] = t1
    B[n-1,n] = 20.0/33.0*hm1
    B[n,n] = 197.0/60.0*hm1

    return A, B
end 

#===============================================
#
#
#
================================================#
function ddT4(n,h)

    dl = fill(0.1,n-1)
    du = fill(0.1,n-1)
    d = ones(n)
    dl[n-1] = 10.0
    du[1] = 10.0
    A = Tridiagonal(dl, d, du)

    B = zeros(n,n)
    hm2 = 1.0/(h*h)
    a = 6.0/5.0
    for i = 2:n-1
        B[i,i-1] = a*hm2
        B[i,i] = -2.0*a*hm2
        B[i,i+1] = a*hm2
    end

    # points 1, N
    B[1,1] = 145.0/12*hm2
    B[1,2] = - 76.0/3*hm2
    B[1,3] = 29.0/2*hm2
    B[1,4] = -4.0/3*hm2
    B[1,5] = 1.0/12*hm2

    B[n,n] = 145.0/12*hm2
    B[n,n-1] = - 76.0/3*hm2
    B[n,n-2] = 29.0/2*hm2
    B[n,n-3] = -4.0/3*hm2
    B[n,n-4] = 1.0/12*hm2

    return A,B
end

#===============================================
#
#
#
================================================#
function ddT6(n,h)
    t1 = 2.0/11.0
    dl = fill(t1,n-1)
    du = fill(t1,n-1)
    dl[1] = 11.0/128.0
    dl[n-2] = 11.0/128.0
    dl[n-1] = 126.0/11.0
    du[1] = 126.0/11.0
    du[2] = 11.0/128.0
    du[n-1] = 11.0/128.0
    d = ones(n)
    A = Tridiagonal(dl, d, du)

    B = zeros(n,n)
    hm2 = 1.0/(h*h)
    t1 = 3.0/44.0
    t2 = 12.0/11.0
    t3 = -2.0*t1 -2.0*t2
    for i = 3:n-2
        B[i,i-2] = t1*hm2
        B[i,i-1] = t2*hm2
        B[i,i] = t3*hm2
        B[i,i+1] = t2*hm2
        B[i,i+2] = t1*hm2
    end

    # points 1, N
    B[1,1] = 2077.0/157*hm2
    B[1,2] = -2943.0/110*hm2
    B[1,3] = 573.0/44*hm2
    B[1,4] = 167.0/99*hm2
    B[1,5] = -18.0/11*hm2
    B[1,6] = 57.0/110*hm2
    B[1,7] = -131.0/1980*hm2

    B[2,1] = 585.0/512*hm2
    B[2,2] = -141.0/64*hm2
    B[2,3] = 459.0/512*hm2
    B[2,4] = 9.0/32*hm2
    B[2,5] = -81.0/512*hm2
    B[2,6] = 3.0/64*hm2
    B[2,7] = -3.0/512*hm2

    B[n-1,n] = 585.0/512*hm2
    B[n-1,n-1] = -141.0/64*hm2
    B[n-1,n-2] = 459.0/512*hm2
    B[n-1,n-3] = 9.0/32*hm2
    B[n-1,n-4] = -81.0/512*hm2
    B[n-1,n-5] = 3.0/64*hm2
    B[n-1,n-6] = -3.0/512*hm2

    B[n,n] = 2077.0/157*hm2
    B[n,n-1] = -2943.0/110*hm2
    B[n,n-2] = 573.0/44*hm2
    B[n,n-3] = 167.0/99*hm2
    B[n,n-4] = -18.0/11*hm2
    B[n,n-5] = 57.0/110*hm2
    B[n,n-6] = -131.0/1980*hm2

    return A,B
end

#===============================================
#
#  These derivative coefficients come from Tables I and II of 
#
#  Jae Wook Kim, "Quasi-disjoint pentadiagonal matrix systems for 
#  the parallelization of compact finite-difference schemes and 
#  filters," Journal of Computational Physics 241 (2013) 168–194.
#
================================================#
function derivKim4(n, h)

    alpha = 0.5862704032801503
    beta = 9.549533555017055e-2

    hm1 = 1.0/h
    a1 = 0.6431406736919156 * hm1
    a2 = 0.2586011023495066 * hm1
    a3 = 7.140953479797375e-3 * hm1

    y00 = 0.0;                   y10=8.360703307833438e-2; 			  y20=3.250008295108466e-2 
    y01=5.912678614078549;       y11= 0.0;				              y21=0.3998040493524358 
    y02=3.775623951744012;       y12=2.058102869495757;               y22=0.0
    y03=0.0;                     y13=0.9704052014790193;              y23=0.7719261277615860 
    y04=0.0;                     y14=0.0;                             y24=0.1626635931256900

    b00=0.0;                       b10=-0.3177447290722621 *hm1;        b20=-0.1219006056449124 *hm1
    b01=-3.456878182643609 *hm1;   b11=0.0;                             b21=-0.6301651351188667 *hm1
    b02=5.839043358834730 *hm1;    b12=-2.807631929593225e-2 *hm1;      b22=0.0
    b03=1.015886726041007 *hm1;    b13=1.593461635747659 *hm1;          b23=0.6521195063966084 *hm1
    b04=-0.2246526470654333 *hm1;  b14=0.2533027046976367 *hm1;         b24=0.3938843551210350 *hm1
    b05=8.564940889936562e-2*hm1;  b15=-3.619652460174756e-2 *hm1;      b25=1.904944407973912e-2 *hm1
    b06=-1.836710059356763e-2*hm1; b16=4.080281419108407e-3 *hm1;       b26=-1.027260523947668e-3*hm1

    b00 = -(b01 + b02 + b03 + b04 + b05 + b06)
    b11 = -(b10 + b12 + b13 + b14 + b15 + b16)
    b22 = -(b20 + b21 + b23 + b24 + b25 + b26)

    P = BandedMatrix(ones(n,n),(2,2))
    P[band(-2)] .= beta
    P[band(-1)] .= alpha
    P[band(1)]  .= alpha
    P[band(2)]  .= beta
    P[1,2] = y01
    P[1,3] = y02
    P[2,1] = y10
    P[2,3] = y12
    P[2,4] = y13
    P[3,1] = y20
    P[3,2] = y21
    P[3,4] = y23
    P[3,5] = y24

    P[n-2,n-4] = y24
    P[n-2,n-3] = y23
    P[n-2,n-1] = y21
    P[n-2,n]   = y20
    P[n-1,n-3] = y13
    P[n-1,n-2] = y12
    P[n-1,n]   = y10
    P[n,n-2]   = y02
    P[n,n-1]   = y01

    Q = BandedMatrix(zeros(n,n),(6,6))
    Q[band(-3)] .= -a3
    Q[band(-2)] .= -a2
    Q[band(-1)] .= -a1
    Q[band(1)]  .= a1
    Q[band(2)]  .= a2
    Q[band(3)]  .= a3

    Q[1,1] = b00
    Q[1,2] = b01
    Q[1,3] = b02
    Q[1,4] = b03
    Q[1,5] = b04
    Q[1,6] = b05
    Q[1,7] = b06

    Q[2,1] = b10
    Q[2,2] = b11
    Q[2,3] = b12
    Q[2,4] = b13
    Q[2,5] = b14
    Q[2,6] = b15
    Q[2,7] = b16

    Q[3,1] = b20
    Q[3,2] = b21
    Q[3,3] = b22
    Q[3,4] = b23
    Q[3,5] = b24
    Q[3,6] = b25
    Q[3,7] = b26

    Q[n-2,n]   = -b20
    Q[n-2,n-1] = -b21
    Q[n-2,n-2] = -b22
    Q[n-2,n-3] = -b23
    Q[n-2,n-4] = -b24
    Q[n-2,n-5] = -b25
    Q[n-2,n-6] = -b26

    Q[n-1,n]   = -b10
    Q[n-1,n-1] = -b11
    Q[n-1,n-2] = -b12
    Q[n-1,n-3] = -b13
    Q[n-1,n-4] = -b14
    Q[n-1,n-5] = -b15
    Q[n-1,n-6] = -b16

    Q[n,n]   = -b00
    Q[n,n-1] = -b01
    Q[n,n-2] = -b02
    Q[n,n-3] = -b03
    Q[n,n-4] = -b04
    Q[n,n-5] = -b05
    Q[n,n-6] = -b06

    return P, Q
end


#===============================================
#
#  These filter coefficients come from Tables III and IV of 
#
#  Jae Wook Kim, "Quasi-disjoint pentadiagonal matrix systems for 
#  the parallelization of compact finite-difference schemes and 
#  filters," Journal of Computational Physics 241 (2013) 168–194.
#
================================================#
function filterKim6(n)

    alphaF=0.6651452077642562
    betaF=0.1669709584471488
    aF1=8.558206326059179e-4
    aF2=-3.423282530423672e-4
    aF3=5.705470884039454e-5
    aF0=-2.0*(aF1+aF2+aF3)

    yF00=0.0;                  yF10=0.7311329755609861;      yF20=0.1681680891936087
    yF01=0.3412746505356879;   yF11=0.0;                     yF21=0.6591595540319565
    yF02=0.2351300295562464;   yF12=0.6689728401317021;      yF22=0.0
    yF03=0.0;                  yF13=0.1959510121583215;      yF23=0.6591595540319565
    yF04=0.0;                  yF14=0.0;                     yF24=0.1681680891936087

    bF20=-2.81516723801634e-4 
    bF21=1.40758361900817e-3
    bF23=2.81516723801634e-3
    bF24=-1.40758361900817e-3
    bF25=2.81516723801634e-4
    bF22= -(bF20 + bF21 + bF23 + bF24 + bF25)

    P = BandedMatrix(ones(n,n),(2,2))
    P[band(-2)] .= betaF
    P[band(-1)] .= alphaF
    P[band(1)]  .= alphaF
    P[band(2)]  .= betaF
    P[1,2] = yF01
    P[1,3] = yF02
    P[2,1] = yF10
    P[2,3] = yF12
    P[2,4] = yF13
    P[3,1] = yF20
    P[3,2] = yF21
    P[3,4] = yF23
    P[3,5] = yF24

    P[n-2,n-4] = yF24
    P[n-2,n-3] = yF23
    P[n-2,n-1] = yF21
    P[n-2,n]   = yF20
    P[n-1,n-3] = yF13
    P[n-1,n-2] = yF12
    P[n-1,n]   = yF10
    P[n,n-2]   = yF02
    P[n,n-1]   = yF01

    Q = BandedMatrix(zeros(n,n),(3,3))
    Q[band(-3)] .= aF3
    Q[band(-2)] .= aF2
    Q[band(-1)] .= aF1
    Q[band(0)]  .= aF0
    Q[band(1)]  .= aF1
    Q[band(2)]  .= aF2
    Q[band(3)]  .= aF3

    Q[1,1] = 0.0
    Q[1,2] = 0.0
    Q[1,3] = 0.0
    Q[1,4] = 0.0

    Q[2,1] = 0.0
    Q[2,2] = 0.0
    Q[2,3] = 0.0
    Q[2,4] = 0.0
    Q[2,5] = 0.0

    Q[3,1] = bF20
    Q[3,2] = bF21
    Q[3,3] = bF22
    Q[3,4] = bF23
    Q[3,5] = bF24
    Q[3,6] = bF25
 
    Q[n-2,n-5] = bF25
    Q[n-2,n-4] = bF24
    Q[n-2,n-3] = bF23
    Q[n-2,n-2] = bF22
    Q[n-2,n-1] = bF21
    Q[n-2,n] = bF20

    Q[n-1,n-4] = 0.0
    Q[n-1,n-3] = 0.0
    Q[n-1,n-2] = 0.0
    Q[n-1,n-1] = 0.0
    Q[n-1,n] = 0.0

    Q[n,n-3] = 0.0
    Q[n,n-2] = 0.0
    Q[n,n-1] = 0.0
    Q[n,n] = 0.0

    return P, Q
end

#===============================================
#
#
#
================================================#
function derivBLT6(n, h)
    hm1 = 1.0/h

    beta_0_p1  = 6.85667616635778

    alpha_0_0  = -3.6546685666048893*hm1
    alpha_0_1  = -2.428065846887595*hm1
    alpha_0_2  = 8.71335233271556*hm1
    alpha_0_3  = -3.5233428330244467*hm1
    alpha_0_4  = 1.03555872211926*hm1

    beta_1_m1  = 0.41230047387468755
    beta_1_p1  = 1.9563380720802086

    alpha_1_0  = -1.0436025117431929*hm1
    alpha_1_1  = 0.0
    alpha_1_2  = -0.7136150600668406*hm1
    alpha_1_3  = 2.3306729849958336*hm1
    alpha_1_4  = -0.6711267770300783*hm1

    beta_2_m1  = -0.47987208402539805
    beta_2_p1  = -1.2930775013841294

    alpha_2_0  = 0.1890770001845506*hm1
    alpha_2_1  = -0.30340795098518447*hm1
    alpha_2_2  = 0.0
    alpha_2_3  = 1.0488462502306883*hm1
    alpha_2_4  = -1.056496112033864*hm1

    beta_3_m1  = -0.3917777288055424
    beta_3_p1  = -0.029222197736104544

    alpha_3_0  = -0.05438332966041568*hm1
    alpha_3_1  = 0.4556295969814727*hm1
    alpha_3_2  = -0.8986296214675904*hm1
    alpha_3_3  = 0.0
    alpha_3_4  = 0.5662870513206056*hm1
    alpha_3_5  = -0.06890369717407233*hm1
    
    gamma_1 = 7.0/9.0*hm1
    gamma_2 = 1.0/36.0*hm1

    delta_1 = 1.0/3.0

    P = BandedMatrix(ones(n,n),(1,1))
    P[band(-1)] .= delta_1
    P[band(1)]  .= delta_1

    P[1,2] = beta_0_p1
    P[2,1] = beta_1_m1
    P[2,3] = beta_1_p1
    P[3,2] = beta_2_m1
    P[3,4] = beta_2_p1
    P[4,3] = beta_3_m1
    P[4,5] = beta_3_p1

    P[n,n-1] = beta_0_p1
    P[n-1,n] = beta_1_m1
    P[n-1,n-2] = beta_1_p1
    P[n-2,n-1] = beta_2_m1
    P[n-2,n-3] = beta_2_p1
    P[n-3,n-2] = beta_3_m1
    P[n-3,n-4] = beta_3_p1

    Q = BandedMatrix(zeros(n,n),(5,5))
    Q[band(-2)] .= -gamma_2
    Q[band(-1)] .= -gamma_1
    Q[band(1)]  .= gamma_1
    Q[band(2)]  .= gamma_2

    Q[1,1] = alpha_0_0
    Q[1,2] = alpha_0_1
    Q[1,3] = alpha_0_2
    Q[1,4] = alpha_0_3
    Q[1,5] = alpha_0_4

    Q[2,1] = alpha_1_0
    Q[2,2] = alpha_1_1
    Q[2,3] = alpha_1_2
    Q[2,4] = alpha_1_3
    Q[2,5] = alpha_1_4

    Q[3,1] = alpha_2_0
    Q[3,2] = alpha_2_1
    Q[3,3] = alpha_2_2
    Q[3,4] = alpha_2_3
    Q[3,5] = alpha_2_4

    Q[4,1] = alpha_3_0
    Q[4,2] = alpha_3_1
    Q[4,3] = alpha_3_2
    Q[4,4] = alpha_3_3
    Q[4,5] = alpha_3_4
    Q[4,6] = alpha_3_5

    Q[n-3,n] = -alpha_3_0
    Q[n-3,n-1] = -alpha_3_1
    Q[n-3,n-2] = -alpha_3_2
    Q[n-3,n-3] = -alpha_3_3
    Q[n-3,n-4] = -alpha_3_4
    Q[n-3,n-5] = -alpha_3_5

    Q[n-2,n] = -alpha_2_0
    Q[n-2,n-1] = -alpha_2_1
    Q[n-2,n-2] = -alpha_2_2
    Q[n-2,n-3] = -alpha_2_3
    Q[n-2,n-4] = -alpha_2_4

    Q[n-1,n] = -alpha_1_0
    Q[n-1,n-1] = -alpha_1_1
    Q[n-1,n-2] = -alpha_1_2
    Q[n-1,n-3] = -alpha_1_3
    Q[n-1,n-4] = -alpha_1_4

    Q[n,n] = -alpha_0_0
    Q[n,n-1] = -alpha_0_1
    Q[n,n-2] = -alpha_0_2
    Q[n,n-3] = -alpha_0_3
    Q[n,n-4] = -alpha_0_4

    return P, Q
end

#===============================================
#
#
#
================================================#
function filterT6(n, alpha)
    dl = fill(alpha, n-1)
    du = fill(alpha, n-1)
    d = ones(n)
    A = Tridiagonal(dl, d, du)
   
    B = zeros(n,n)
   
    a = (11 + 10*alpha)/16
    b = (15 + 34*alpha)/32
    c = (-3 + 6*alpha)/16
    d = (1 - 2*alpha)/32

    for i = 4:n-3
        B[i,i-3] = d/2.0
        B[i,i-2] = c/2.0
        B[i,i-1] = b/2.0
        B[i,i] = a
        B[i,i+1] = b/2.0
        B[i,i+2] = c/2.0
        B[i,i+3] = d/2.0
    end
   
    B[1,1] = (63 + alpha) / 64
    B[1,2] = (3 + 29 * alpha) / 32
    B[1,3] = 15 * (-1 + alpha) / 64
    B[1,4] = 5 * (1 - alpha) / 16
    B[1,5] = 15 * (-1 + alpha) / 64
    B[1,6] = 3 * (1 - alpha) / 32
    B[1,7] = (-1 + alpha) / 64
   
    B[2,1] = (1 + 62 * alpha) / 64
    B[2,2] = (29 + 6 * alpha) / 32
    B[2,3] = (15 + 34 * alpha) / 64
    B[2,4] = (-5 + 10 * alpha) / 16
    B[2,5] = 15 * (1 - 2 * alpha) / 64
    B[2,6] = 3 * (-1 + 2 * alpha) / 32
    B[2,7] = (1 - 2 * alpha) / 64

    B[3,1] = (-1 + 2 * alpha) / 64
    B[3,2] = (3 + 26 * alpha) / 32
    B[3,3] = (49 + 30 * alpha) / 64
    B[3,4] = (5 + 6 * alpha) / 16
    B[3,5] = 15 * (-1 + 2 * alpha) / 64
    B[3,6] = 3 * (1 - 2 * alpha) / 32
    B[3,7] = (-1 + 2 * alpha) / 64

    B[n-2,n] = (-1 + 2 * alpha) / 64
    B[n-2,n-1] = (3 + 26 * alpha) / 32
    B[n-2,n-2] = (49 + 30 * alpha) / 64
    B[n-2,n-3] = (5 + 6 * alpha) / 16
    B[n-2,n-4] = 15 * (-1 + 2 * alpha) / 64
    B[n-2,n-5] = 3 * (1 - 2 * alpha) / 32
    B[n-2,n-6] = (-1 + 2 * alpha) / 64

    B[n-1,n] = (1 + 62 * alpha) / 64
    B[n-1,n-1] = (29 + 6 * alpha) / 32
    B[n-1,n-2] = (15 + 34 * alpha) / 64
    B[n-1,n-3] = (-5 + 10 * alpha) / 16
    B[n-1,n-4] = 15 * (1 - 2 * alpha) / 64
    B[n-1,n-5] = 3 * (-1 + 2 * alpha) / 32
    B[n-1,n-6] = (1 - 2 * alpha) / 64

    B[n,n] = (63 + alpha) / 64
    B[n,n-1] = (3 + 29 * alpha) / 32
    B[n,n-2] = 15 * (-1 + alpha) / 64
    B[n,n-3] = 5 * (1 - alpha) / 16
    B[n,n-4] = 15 * (-1 + alpha) / 64
    B[n,n-5] = 3 * (1 - alpha) / 32
    B[n,n-6] = (-1 + alpha) / 64

    return A, B
end

#===============================================
#
#
#
================================================#
function filterT8(n, alpha)
    dl = fill(alpha, n-1)
    du = fill(alpha, n-1)
    d = ones(n)
    A = Tridiagonal(dl, d, du)

    B = zeros(n,n)

    a = (93+70*alpha)/128
    b = (7+18*alpha)/16
    c = (-7 + 14*alpha)/32
    d = (1 -2*alpha)/16
    e = (-1 + 2*alpha)/128

    for i = 5:n-4
        B[i,i-4] = e/2.0
        B[i,i-3] = d/2.0
        B[i,i-2] = c/2.0
        B[i,i-1] = b/2.0
        B[i,i] = a
        B[i,i+1] = b/2.0
        B[i,i+2] = c/2.0
        B[i,i+3] = d/2.0
        B[i,i+4] = e/2.0
    end

    B[1,1] = (255 + alpha) / 256
    B[1,2] = (1 + 31 * alpha) / 32
    B[1,3] = (-7 + 7 * alpha) / 64
    B[1,4] = (7 - 7 * alpha) / 32
    B[1,5] = 7 * (-5 + 5 * alpha) / 128
    B[1,6] = (7 - 7 * alpha) / 32
    B[1,7] = 7 * (-1 + alpha) / 64
    B[1,8] = (1 - alpha) / 32
    B[1,9] = (-1 + alpha) / 256

    B[2,1] = (1 + 254 * alpha) / 256
    B[2,2] = (31 + 2 * alpha) / 32
    B[2,3] = (7 + 50 * alpha) / 64
    B[2,4] = (-7 + 14 * alpha) / 32
    B[2,5] = 7 * (5 - 10 * alpha) / 128
    B[2,6] = (-7 + 14 * alpha) / 32
    B[2,7] = (7 - 14 * alpha) / 64
    B[2,8] = (-1 + 2 * alpha) / 32
    B[2,9] = (1 - 2 * alpha) / 256

    B[3,1] = (-1 + 2 * alpha) / 256
    B[3,2] = (1 + 30 * alpha) / 32
    B[3,3] = (57 + 14 * alpha) / 64
    B[3,4] = (7 + 18 * alpha) / 32
    B[3,5] = 7 * (-5 + 10 * alpha) / 128
    B[3,6] = (7 - 14 * alpha) / 32
    B[3,7] = (-7 + 14 * alpha) / 64
    B[3,8] = (1 - 2 * alpha) / 32
    B[3,9] = (-1 + 2 * alpha) / 256

    B[4,1] = (1 - 2 * alpha) / 256
    B[4,2] = (-1 + 2 * alpha) / 32
    B[4,3] = (7 + 50 * alpha) / 64
    B[4,4] = (25 + 14 * alpha) / 32
    B[4,5] = (35 + 58 * alpha) / 128
    B[4,6] = (-7 + 14 * alpha) / 32
    B[4,7] = (7 - 14 * alpha) / 64
    B[4,8] = (-1 + 2 * alpha) / 32
    B[4,9] = (1 - 2 * alpha) / 256

    B[n-3,n] = (1 - 2 * alpha) / 256
    B[n-3,n-1] = (-1 + 2 * alpha) / 32
    B[n-3,n-2] = (7 + 50 * alpha) / 64
    B[n-3,n-3] = (25 + 14 * alpha) / 32
    B[n-3,n-4] = (35 + 58 * alpha) / 128
    B[n-3,n-5] = (-7 + 14 * alpha) / 32
    B[n-3,n-6] = (7 - 14 * alpha) / 64
    B[n-3,n-7] = (-1 + 2 * alpha) / 32
    B[n-3,n-8] = (1 - 2 * alpha) / 256

    B[n-2,n] = (-1 + 2 * alpha) / 256
    B[n-2,n-1] = (1 + 30 * alpha) / 32
    B[n-2,n-2] = (57 + 14 * alpha) / 64
    B[n-2,n-3] = (7 + 18 * alpha) / 32
    B[n-2,n-4] = 7 * (-5 + 10 * alpha) / 128
    B[n-2,n-5] = (7 - 14 * alpha) / 32
    B[n-2,n-6] = (-7 + 14 * alpha) / 64
    B[n-2,n-7] = (1 - 2 * alpha) / 32
    B[n-2,n-8] = (-1 + 2 * alpha) / 256

    B[n-1,n] = (1 + 254 * alpha) / 256
    B[n-1,n-1] = (31 + 2 * alpha) / 32
    B[n-1,n-2] = (7 + 50 * alpha) / 64
    B[n-1,n-3] = (-7 + 14 * alpha) / 32
    B[n-1,n-4] = 7 * (5 - 10 * alpha) / 128
    B[n-1,n-5] = (-7 + 14 * alpha) / 32
    B[n-1,n-6] = (7 - 14 * alpha) / 64
    B[n-1,n-7] = (-1 + 2 * alpha) / 32
    B[n-1,n-8] = (1 - 2 * alpha) / 256

    B[n,n] = (255 + alpha) / 256
    B[n,n-1] = (1 + 31 * alpha) / 32
    B[n,n-2] = (-7 + 7 * alpha) / 64
    B[n,n-3] = (7 - 7 * alpha) / 32
    B[n,n-4] = 7 * (-5 + 5 * alpha) / 128
    B[n,n-5] = (7 - 7 * alpha) / 32
    B[n,n-6] = 7 * (-1 + alpha) / 64
    B[n,n-7] = (1 - alpha) / 32
    B[n,n-8] = (-1 + alpha) / 256

    return A, B
end


#===============================================
#
#
#
================================================#
function debugL6(dxu, u, h)

    s = 1
    p = 2



end

#===============================================
#
#
#
================================================#
