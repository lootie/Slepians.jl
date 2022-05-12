using Slepians, FastGaussQuadrature, Dierckx

"""

    scale_quad_nodez(Nqz, z)
    
Scale the quadrature nodes in z to the correct interval

# Arguments
- `Nqz`: Number of quadrature nodes in z
- `z`: Levels in z given

# Outputs
- `colf`: Zero crossing indices where the thph z-coordinate and quadrature node intersect
- `th`: z-level at which crossings are required
- `wqz`: Quadrature weights in the z-direction

"""
function scale_quad_nodez(Nqz, z)
    # Calculate blanks, z-axis GL points regardless of the data range
    qz, wqz = FastGaussQuadrature.gausslegendre(Nqz)

    mn, mx = (minimum(z), maximum(z))
    # Scale the z-quadrature points to the min-max interval
    th = qz*(abs(mx - mn)/2) .+ (mx + mn)/2;

    # We need to find the two z's between which each quadrature point (th) lives
    xmth = repeat(z, 1, length(th)) .- repeat(th', length(z), 1)
    dsx = diff(sign.(xmth), dims = 1)
    col = findall(x -> (x == 2), dsx)
    colf = map(c -> Tuple(c)[1], col)
    return colf, th, wqz
end

"""

    fill2d(colf, thpha, th, Nqx, Nqy, N)
    
Fill each 2D level with nodes in x and y

# Arguments
- `colf`: Zero crossing indices where the thph z-coordinate and quadrature node intersect
- `thpha`: The array of arrays of closed curves where the columns are in order (z, y, x)
- `th`: z-level at which crossings are required
- `Nqx`: Number of quadrature nodes in x
- `Nqy`: Number of quadrature nodes in y
- `N::Int64`: Number of points in each closed curve

# Outputs
- `pkg`: The array containing quadrature nodes in x, y at each level z

"""
function fill2d(z, colf, thpha, th, Nqx, Nqy)
    N = size(thpha[1], 1)
    newcurve = map(i -> interpcontour(Float64.(z[colf[i]:(colf[i] + 1)]), th[i], vcat(thpha[colf[i]], thpha[(colf[i] + 1)]), N), 1:length(th))
    pkg = map(nc -> get_quadrature_nodes_2D(nc[:,3], nc[:,2], Nqx, Nqy), newcurve);
    return newcurve, pkg
end

"""

    getnodeswts3d(szs, pkg, th, wqz)
    
Get the nodes and weights in 3D

# Arguments
- `szs::Tuple{Int64}`: number of quadrature nodes in each dimension
- `pkg`: The array containing quadrature nodes in x, y at each level z
- `th`: z-level at which crossings are required 
- `wqz`: Quadrature weights in the z-direction

# Outputs
- `no`: The quadrature nodes in Tuple format
- `sqwt`: squared weight associated with each node
- `ev`: The even grid of size szs

"""
function getnodeswts3d(szs, pkg, th, wqz)
    QX = mapreduce(p -> p[1][:], vcat, pkg)
    QY = mapreduce(p -> p[2][:], vcat, pkg)
    QZ = mapreduce(i -> ones(size(pkg[i][2][:]))*th[i], vcat, 1:length(pkg)); #kron(th, ones(Nqx, Nqy))[:]

    wts = mapreduce(i -> ones(size(pkg[i][2][:]))*wqz[i], vcat, 1:length(pkg));

    no    = map(i->(QX[i], QY[i], QZ[i]), 1:prod(size(QX)))
    sqwt = sqrt.(wts);
    ev = vec(collect(product([range(minimum(QX), maximum(QX), length = szs[1]), range(minimum(QY), maximum(QY), length = szs[2]), range(minimum(QZ), maximum(QZ), length = szs[3]) ]...)))

    return no, sqwt, ev
end

""" 
    respline(x, y, N)

Use a 2D parametric B-spline to interpolate the closed curve to N points

# Arguments
- `x`: vector of x-coordinates
- `y`: vector of y-coordinates
- `N`: number of desired output points

# Outputs
- a matrix of size 2 x N in [y, x] order containing the splined coordinates

"""
function respline(x, y, N)
    ps = ParametricSpline(vcat(y', x'); s=0.0, k=1)
    return evaluate(ps,  LinRange(0, 1, N))
end

""" Spline all curves to have the same number of points """
function equalNclosedcurve(thpha, N) 
    output = []
    for thph in thpha
        if (length(thph) == N)
            push!(output, thph) 
        else
            try push!(output, hcat(thph[1,1]*ones(N), respline(thph[:,3], thph[:,2], N)')) catch; [] end
            #push!(output, hcat(thph[1,1]*ones(N), res))
        end
    end
    return output
end

"""

    slepian3(szs, N, z, thpha; <kwargs>)
    
Compute the 3D localization problem with boundary defined as a set of stacked contours

# Arguments
- `szs::Tuple{Int64}`: number of quadrature nodes in each dimension
- `z`: Levels in z given. These must correspond to the order in thpha
- `thpha`: The array of arrays of closed curves where the columns are in order (z, y, x)

# Optional Keyword Arguments
- `M::Int64 = 3`: Number of Slepians to output
- `Kp::Vector{Float64}`: Bandwidth (radius) in units of number of pixels
- `exact::Bool`: Whether or not to use the approximation using HODLR matrices (false) or not (true)
- `prec::Float63`: If exact = false, state precision of the approximation
- `lvl::Int64`: Number of levels of the HODLR approximation
- `maxrank::Int64`: Maximum rank of the off-diagonal matrices
- `int::Tuple{Int64}`: output size (if different than szs) to interpolate to (Nyström)

# Outputs
- `s`: 3D Slepians
- `sl`: eigenvalues

"""
function slepian3(szs, z, thpha; M = 3, Kp = [4.0], exact = false,
                    prec = 1e-8, lvl = 6, 
                    maxrank = 256, int = szs)
    Nqx, Nqy, Nqz = szs 
    # compute the nodes - first scale in z
    colf, th, wqz = scale_quad_nodez(Nqz, z)
    # Make sure the closed contours have the same number of points
    N = maximum(size.(thpha,1))
    thpha = equalNclosedcurve(thpha, N)
    
    # compute the quad nodes in 2d at each level
    pkg = fill2d(z, colf, thpha, th, Nqx, Nqy)[2]
    
    # collate nodes weights
    no, sqwt, ev = getnodeswts3d(szs, pkg, th, wqz)
    
    # do the integration
    s, sl = customsleps_ext(M, Kp, szs, prec = prec, exact = exact, lvl = lvl, maxrank = maxrank, no = no,
            sqwt = sqwt, int = nothing, ev = ev);
    return s, sl
end
        
# Compute the magnitude squared FFT of each component of the vector of matrices sl
SL(sl) = map(i->abs2.(fftshift(fft(sl[i]))), 1:length(sl))   