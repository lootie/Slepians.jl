"""
    getbdypts_2d(mask)

Starting with a matrix of zeros and ones, get the points on the boundary

# Arguments
    - `mask` a matrix of zeros and ones

# Outputs
    - `x` x-coordinates of the boundary points
    - `y` y-coordinates of the boundary points
"""
function getbdypts_2d(mask::Array{T,2}) where T<:Number
    # difference rows and columns
    dif1 = diff(mask, dims = 1)
    dif2 = diff(mask, dims = 2)
    # find the nonzero entries of the row differences (up, down respectively)
    col1u = findall(dif1 .> 0)
    col1d = findall(dif1 .< 0)

    T1, T2 = size(mask)
    t1 = 1:T1
    t2 = 1:T2

    # this translates into (-1, 1) distance from index
    x_up = repeat(t1', length(t2), 1)
    y_up = repeat(t2, 1, length(t1))

    # find the nonzero entries of the column differences (up, down respectively)
    col2u = findall(dif2 .> 0)
    col2d = findall(dif2 .< 0)

    xp = vcat(x_up[col1u], x_up[col1d], x_up[col2u], x_up[col2d])
    yp = vcat(y_up[col1u], y_up[col1d], y_up[col2u], y_up[col2d])

    per = sortperm(xp)

    return xp[per], yp[per]
end


""" 
    findnst(xs, ys, x0, y0)

Find the nearest 2D point in arrays xs, ys to (x0, y0) Euclidean distance 
"""
function findnst(xs::Vector, ys::Vector, x0, y0)
    ind = findmin((xs .- x0).^2 .+ (ys .- y0).^2)[2]
    return (ind, xs[ind], ys[ind])
end

"""
    closedcurve_2d(xs, ys)

Generate a closed curve based on (x, y) points given and Euclidean distance

# Arguments
- `xs` x coordinates of the points
- `ys` y coordinates of the points

# Outputs
- `xpath` x coordinates of the path
- `ypath` y coordinates of the path
"""
function closedcurve_2d(xs, ys)
    ind = zeros(Int64, length(xs))
    xpath = zeros(length(xs))
    xpath[1] = xs[1]
    ypath = zeros(length(ys))
    ypath[1] = ys[1]
    compset = collect(1:length(xs))
    for i = 2:length(xs)
        ind[i], xpath[i], ypath[i] = findnst(xs[compset], ys[compset], xpath[i-1], ypath[i-1])
        # xpath[i], ypath[i] = (xs[ind[i]], ys[ind[i]])
        compset = compset[1:end .!= ind[i]]
    end
    push!(xpath, xpath[1])
    push!(ypath, ypath[1])
    return xpath, ypath
end

function mask2closedcurve(mask)
    xs, ys = getbdypts_2d(Float64.(mask))
    return closedcurve_2d(xs, ys)
end
