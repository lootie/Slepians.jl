using Slepians, FastGaussQuadrature
using Plots, FFTW

M = 4; Kp = 10.0; szs = (16,16,16)
int = (51, 51, 51)
nowtv = [FastGaussQuadrature.gausslegendre(nj) for nj in szs]

lbfun(x) = (x >= 0) * max(-sqrt(2) + x, -1) + (x < 0) * max(-sqrt(2) - x, -1)

no    = [(x, lbfun(x)*y, z) for x in nowtv[1][1] for y in nowtv[2][1] for z in nowtv[3][1]];

using Profile

# scatter(no, camera = (60, 60))

# Compute the Slepians on this basis
lvl = 6
maxrank = 128
@time s, sleps = customsleps(M, Kp, szs, no = no, sqwt = givewts(nowtv), int = int, 
    exact = false, lvl = lvl, maxrank = maxrank);

using HDF5

A = collect(reshape(1:120, 15, 8))

h5open("slep100.h5", "w") do file
    write(file, "s", s)  
    write(file, "M", M)  
    write(file, "Kp", Kp)
    write(file, "szs", szs)
    for i in 1:M
        write(file, "slep$i", sleps[i])
    end
end

#=
# To read in the file:
filename = "slep100.h5"
using HDF5 

M = h5open(filename, "r") do file
    read(file, "M")
end

s = h5open(filename, "r") do file
    read(file, "s")
end

sleps = h5open(filename, "r") do file
    map(i -> read(file, "slep$i"), 1:M)
end
=#


