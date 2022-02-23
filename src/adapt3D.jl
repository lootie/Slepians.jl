""" 
    interp2(z, xy, z0)

# Arguments
- `z::Vector` vector of z's
- `yx::Matrix` y, x coordinates at each of the z's
- `z0<:Number` level at which to interpolate z0

# Outputs
- zyx coordinates corresponding to level z0

"""
interp2(z::Vector, yx::Matrix, z0 = 0.0) = (yx[2,:] - yx[1,:])*(z0 - z[1])/(z[2] - z[1]) .+ yx[1,:]


# To get the closed curve at the level th, take the 100 points before and 100 points after 
# the crossing and interpolate each yx
"""
    interpcontour(z, z0, thph, N)

# Arguments
- `z::Vector` the two z-values between which the level z0 falls. 
- `z0<:Number` the z-level of the desired contour
- `thph::Matrix` the zyx coordinates (first N are at level z[1] second are at level z[2])
- `N::Int64` number of points for each contour
"""
function interpcontour(z, z0, thph, N)
    newcurve = map(i -> interp2(z, thph[[i, i + N],:], z0), 1:N) 
    return hcat(newcurve...)'
end
