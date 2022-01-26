
# Software note: this code was originally developed by F. Simons and D. Wang in
# the matlab programming language, and was rewritten in julia for performance.
# The original code is available at:
# https://github.com/csdms-contrib/slepian_alpha
# and is distributed under the GNU GPL v 2.0 license.
# I take responsibility for any errors in this code as it is not a verbatim copy
# of the original. In particular, some inbuilt matlab functions (e.g. sub2ind,
# interp1) have been replaced with Julia ones. 

# This software is licensed to C. L. Haley in 2021 under GNU GPL v 2.0

using Interpolations, PCHIPInterpolation

""" 

    sub2ind(A, row, col)

Convert to linear indices 
https://www.mathworks.com/help/matlab/ref/sub2ind.html
"""
function sub2ind(A, row, col) 
    LinearIndices(A)[CartesianIndex(row, col)]
end



"""

    interp1()

1D data interpolation, linear
https://www.mathworks.com/help/matlab/ref/interp1.html?s_tid=doc_ta

# Arguments
- Vector x contains the sample points, and 
- v contains the corresponding values, v(x). 
- Vector xq contains the coordinates of the query points.

# Outputs
- interpolated values of a 1-D function at specific query points using linear interpolation

# Example usage

```
interp1(LinRange(0, 1, 10), log10.(LinRange(0, 1, 10)), [0.5, 0.9], :linear)
```

"""
function interp1(x, v, xq, method = :linear)
    #(method != :linear)||(method != :pchip) && error("Linear Interpolation is the only one implemented here.")
    interpolant = (method == :linear) ? LinearInterpolation(x, v) : PCHIPInterpolation.Interpolator(x, v)
    return interpolant.(xq)
end

"""

    gamini(data, folding)

# Arguments
-`data`: some data vector
-`folding`: The replication factor for every element of the data; if a scalar thsi applies to all of the elements (default 3), if zero or negative, no replication occurs

# Outputs
-`bigger`

# Example usage
```
a, b = degamini(gamini([1, 2, 3, 1, 4, 5], [1, 2, 3, 2, 4, 2]))
```
One gets [1, 2, 2, 3, 3, 3, 1, 1, 4, 4, 4, 4, 5, 5] as the intermediate result.

# See also
@degamini, @gamini2
"""
function gamini(data, folding = 3)
    # Copy a single scalar folding for all elements in the array
    if prod(size(folding)) == 1
        folding = repeat(folding, prod(size(data)), 1)
    end
    # Never had a replication factor of zero before
    # But only if they are of the same length!
    if length(data) == length(folding)
        ind = findall(f->f>0, folding)
        data = data[ind]
        folding = folding[ind]
    end
    # Rearrange
    data = data[:]
    folding = folding[:]
    !(size(data) == size(folding)) && error("Sizes of input and folding must be the same.")
    gelp = zeros(Int64, sum(folding))
    gelp[vcat(1, cumsum(folding[1:end-1], dims = 1) .+ 1)] .= 1
    if !isempty(data)
        bigger  = data[cumsum(gelp, dims = 1)]
    else
        bigger = []
    end
    return bigger
end

""" 

    degamini(v) 

# Arguments
- `v` A vector with repeated entries

# Outputs
- `dv` The same vector with all the repeat sset to 1
- `foldi` A same-dimensional vector with how many repeats there are
- `be` A matrix with begin and end indices into the original vector

# Example
```
dv, foldi, be = degamini([1, 2, 2, 3])
```
answer is dv = [1, 2, 3]; foldi = [1, 2, 1]; be = [1 1; 2 3; 4 4]
"""
function degamini(v)
    v = v[:]
    indi = vcat(1, findall(diff(v) .!= 0) .+ 1) 
    dv = v[indi]
    foldi = diff(vcat(indi, length(v) + 1))
    (length(v) != sum(foldi)) && error("Incorrect")
    be = hcat(indi[:], cumsum(foldi[:]))
    return dv, foldi, be
end

""" 
    
    matranges(ranges)
    
Makes an index vector with monotonically increasing indices between pairs of numbers supplied as input.
From slepian_alpha
    
# Arguments
- `ranges`

    # Outputs
    
# Example
```
matranges([1, 4, 1, 2, -1, 2])
# answer is [1, 2, 3, 4, 1, 2, -1, 0, 1, 2]
```
"""
function matranges(ranges)
    ranges = ranges[:]
    (mod(length(ranges), 2) != 0) && error("Ranges must form pairs")
    lower = ranges[1:2:end]
    upper = ranges[2:2:end]
    hulp1 = ones(sum(upper .- lower .+ 1))
    hulp2 = cumsum(upper[1:end-1] .- lower[1:end-1] .+ 1)
    hulp1[hulp2 .+ 1] = lower[2:end] .- upper[1:(end-1)]   
    hulp1[1] = ranges[1]
    return Int64.(cumsum(hulp1))
end

""" 

    randcirc(xm, ym, r,dr, N)

# Arguments
- xm horizontal positon of the center
- ym vertical position of the center
- `r` radius
- `dr` size of random perturbations around the radius
- `N` number of random spike points

# Outputs
- `x` x-coordinate
- `y` y-coordinate

# Related
@blob

"""
function randcirc(xm = 0.0, ym = 0.0, r = 1.0, dr = 1, N = 10)
    nr = 100
    r = r*ones(N)
    r = r + 2dr*(rand(N) .- 0.5)
    t = LinRange(0, 2*pi*N/(N+1), N)
    
    r = vcat(r,r,r)
    t = vcat(t .- 2 * pi, t, t .+ 2 *pi)
    
    tt = LinRange(0, 2*pi, nr)
    rr = interp1(t, r, collect(tt), :pchip)
    
    x = @. xm + rr * cos(tt)
    y = @. ym + rr * sin(tt)
    
    return x, y
end

""" 
    blob(N, Nj)

Makes moving picture of a random blob by superposition of random circles

# Arguments
- `N` number of loops for, and if movie, default = 100
- `Nj` smoothness, roughly (?)

# Outputs
- `x`
- `y` 
"""
function blob(N = 100, Nj = 10)
    xold, yold = randcirc(0, 0, 1, 1, 10) ### Help: what is randcirc
    r = LinRange(0, 1, Nj+1)
    rm = 1 .- r
    xx = []
    yy = []
    for index = 1:N
        x,y = randcirc(0, 0, 1, 0.2, 10)
        for j = 1:Nj
            push!(xx, x * r[j] + xold * rm[j])
            push!(yy, y * r[j] + yold * rm[j])
        end
        xold = x
        yold = y
    end
    return xx, yy
end

""" 

    phicurve(thph, th)

Adapted from slepian_alpha; finds the longitude crossings and thus integration domains of a closed curve parameterized in colatitude/longitude space at certian query points of colatitude.

# Arguments
- `thph::` Colatitude/longitude of the closed curve (degrees)
- `th::` Colatitude at which crossings are required (degrees)

# Outputs
- `phint` A matrix with crossings/intervals and zeros of dimensions MxN where M = length(th) and N can be anything depending on the oscillations of the curve
- `thp` Colatitude matrix for hatched plotting, if possible
- `php` Longitude matrix for hatched plotting, if possible
- `forreal` Indices of the ones that are real (could be at zero)

Depends on: @sub2ind, @interp1, @degamini, @matranges, and possibly @blob (demo2)
"""
function phicurve(thph, th)
    # Fore every th, find the relevant phint
    xmth = repeat(thph[:,1], 1, length(th)) .- repeat(th', length(thph[:,1]), 1)
    dsx = diff(sign.(xmth), dims = 1)
    # It can be the one before, or after the crossing
    # colf, colj = findall(x -> (x != 0.0), dsx)
    col = findall(x -> (x != 0.0), dsx)
    colf, colj = (map(c -> Tuple(c)[1], col), map(c -> Tuple(c)[2], col))
    # colr = sub2ind(dsx, colf, colj) 
    colr = LinearIndices(dsx)[col]
    # This returns the one on the negative side of hte line
    # add one to colx if the difference is -2; add one to colx2 if the difference is 2
    colx = colf .+ (dsx[colr] .== -2) # DOT was missing here in parentheses
    colx2 = colf .+ (dsx[colr] .== 2)
    L = length(colx)
    (L%2 == 1) && error("Cannot find pairs of crossings.")
    phint = zeros(L)
    # Then one point was exactly hit, this is the thN or the thS case
    if length(colx) == 2 && colx == colx2
        phint = thph[hcat(colx[2], colx2[2]), 2]
        thp = hcat(th, th)
        php = phint
    else
        for ond = 1:L
            # In case you have a node immediately followed by a crossing
            phint[ond] = (colx[ond] == colx2[ond]) ? NaN : interp1(xmth[hcat(colx[ond], colx2[ond]), colj[ond]][:], thph[hcat(colx[ond], colx2[ond]),2][:], 0, :linear) 
        end
    end
    # If the NaNs are not consecutive pairs, get special case
    # Now rearrange back othte number of requested points
    # There could be points with more or less than two crossings
    # Maximum number of times a crossing is repeated
    a, b =  degamini(colj) 
    rowj = copy(colj)
    colj = matranges(Int64.(reshape(hcat(ones(length(b)), b)', length(b)*2, 1)))
    pint = repeat([NaN], length(th), maximum(b))
    subsi = (colj .- 1) .* length(th) .+ rowj # Linear index
    pint[subsi] = phint
    wt, thp = (length(b) == length(th)) ? (0, reshape(gamini(th, b), (2, Int64(length(phint)/2)))) : (1, [])
    # Need to sort since contour may be given in any order
    phint = sort(pint, dims = 2)
    php = (wt == 0) ? reshape(phint[subsi], (2, Int64(length(colj)/2))) : [] # line 95
    # Make them zero so the integral doesn't do anything
    forreal = map(x -> !isnan(x), phint)
    phint[findall(isnan, phint)] .= 0.0
    # note: can use (demo2)
    # x,y=blob(1,1); thph=hcat(y[:],x[:]); Nth = ceil(rand*300); th=LinRange(minimum(thph[:,1]), maximum(thph[:,1]), Nth)
    return phint, thp, php, forreal
end

""" Scale the quaduature points and weights """
function quadpts(qx, Nqx, qy, forreal, xints, wqx, wqy)
    Nall = Int64(sum(forreal[:])/2)
    # Initialize quadrature points
    QX = repeat([NaN], Nqx, Nall)
    QY = repeat([NaN], Nqx, Nall)
    w = repeat([NaN], Nqx, Nall)
    Nrun = 0
    for yindex = 1:length(qy)
        # How many horizontal intervals at this y?
        Nxint = Int64(sum(forreal[yindex,:])/2)
        abset = reshape(xints[yindex, 1:Nxint*2], 2, Nxint)
        # beginning and endpoints of those intervals
        a = abset[1,:]
        b = abset[2,:]
        # produce the scaled nodes, a column per interval
        exs = repeat(a', length(qx), 1) .+ (qx .+ 1)/2 .* (b .- a)'
        # scatter!(p, exs, -0.5*ones(32), marker = :o)
        # And produce the weights that go with it 
        # and multiply it with the current y node and weight 
        # and produce copies of the y nodes, as many as needed 
        # and put them into two big vectors of weights and points
        QX[:, (Nrun+1):(Nrun+Nxint)]= exs
        QY[:, (Nrun+1):(Nrun+Nxint)] = repeat([qy[yindex]], size(exs)...)
        w[:, (Nrun+1):(Nrun+Nxint)] = wqx.*((b.-a)'/2)*wqy[yindex]
        Nrun = Nrun + Nxint
    end
    return QX, QY, w, Nrun
end

sign(x) = (x >= 0.0) ? 1 : -1

function get_quadrature_nodes_2D(x, y, Nqx = 32, Nqy = 32)
    #println("Gauss-Legendre method with ($Nqx, $Nqy) integration nodes.")

    # Calculate blanks, xaxis GL points regardless of the data range
    qx, wqx = FastGaussQuadrature.gausslegendre(Nqx)
    qy, wqy = FastGaussQuadrature.gausslegendre(Nqy)
    
    mn, mx = (minimum(y), maximum(y))
    # Scale the y-quadrature points to the min-max interval
    th = qy*(abs(mx - mn)/2) .+ (mx + mn)/2

    # For each of these points in y find the appropriate ranges of x's
    xints,yp,xp,forreal = phicurve(hcat(y, x), th) 
    QX, QY, w, Nrun = quadpts(qx, Nqx, th, forreal, xints, wqx, wqy)

    return QX, QY, w, Nrun
end