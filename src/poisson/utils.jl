#-----------------------------------------------------------------------
# This function defines coordinate arrays.  
# The coordinate x ranges from xmin to xmax with n points.
#-----------------------------------------------------------------------
function defineCoords(n, xmin, xmax)
    x = Array{Float64}(undef,n)
    dx = (xmax - xmin)/(n-1)
    for i = 1:n
        x[i] = xmin + (i-1)*dx
    end
    return x
end

#-----------------------------------------------------------------------
# Take x-derivative of a function with 2nd order precision, 1st on boundaries.
#-----------------------------------------------------------------------
function deriv_x21(df, f, h)
    nx, ny = size(f)
    hm1 = 1.0/h
    for j = 1:ny
        df[1,j] = (f[2,j] - f[1,j])*hm1
        for i = 2:nx-1
            df[i,j] = 0.5*(f[i+1,j] - f[i-1,j])*hm1
        end
        df[nx,j] = (f[nx,j] - f[nx-1,j])*hm1
    end
end

#-----------------------------------------------------------------------
# Take y-derivative of a function with 2nd order precision, 1st on boundaries.
#-----------------------------------------------------------------------
function deriv_y21(df, f, h)
    nx, ny = size(f)
    hm1 = 1.0/h
    for i = 1:nx
        df[i,1] = (f[i,2] - f[i,1])*hm1
        for j = 2:ny-1
            df[i,j] = 0.5*(f[i,j+1] - f[i,j-1])*hm1
        end
        df[i,ny] = (f[i,ny] - f[i,ny-1])*hm1
    end
end

#-----------------------------------------------------------------------
# Thin a 2D array.  Used for quiver plots.  Takes every nth point in
# each direction from array.
#-----------------------------------------------------------------------
function thindata(f, nth)
    nx0, ny0 = size(f)
    nx = trunc(Int, (nx0-1)/nth) + 1
    ny = trunc(Int, (ny0-1)/nth) + 1
    fn = zeros(nx, ny)
    for j=1:ny
        for i = 1:nx
            fn[i,j] = f[(i-1)*nth + 1,(j-1)*nth + 1]
        end
    end
    return fn
end


#-----------------------------------------------------------------------
# Calculate the norm of 2d vector function
#-----------------------------------------------------------------------
function vectorFieldNorm(vx, vy)
    Enorm = similar(vx)
    @. Enorm = sqrt(vx^2 + vy^2) + 1.0e-8
end


