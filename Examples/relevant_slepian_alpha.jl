
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
function interp1(x, v, xq, method)
    (method != :linear)||(method != :pchip) && error("Linear Interpolation is the only one implemented here.")
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
function randcirc(xm = 0.0, ym = 0.0, r = 1.0, dr = 01, N = 10)
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

Makes moving picture of a random blob by superposition of random sircles

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
