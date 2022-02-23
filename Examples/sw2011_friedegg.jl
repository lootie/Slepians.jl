using MATLAB, Plots

# Make sure to have slepian_alpha and slepian_foxtrot repos available to matlab.
# I did this by exporting the path in my bashrc as follows:
# export MATLABPATH=$HOME/Repos/slepian_alpha:$HOME/Repos/slepian_foxtrot
# you will need to have my forked versions of the repos (lootie fork)

# This recreates figure four of the Simons and Wang Paper. 

method = mxarray("SE")
circn = 256
XY = hcat(cos.(LinRange(0, 2*pi, circn)), sin.(LinRange(0, 2*pi, circn)))

# mdis = mean(sum(sqrt(diff(XY).^2), dims = 2))

N = 42
J = 30
j = 7
xp = yp = LinRange(-2, 2.0, 2^j)

Nqx = Nqy = 32
dxdy = []
XYP = hcat(repeat(xp, 1, length(yp))[:], repeat(yp', length(xp))[:])

G, H, V, K, XYP, XY, A = mxcall(:swfried_easy, 7, mxarray("GL"))
# G, H, V, K, XYP, XY, A, c11, cmn = mxcall(:localization2D, 9, XY, N, J, [], Nqx, Nqy, dxdy, XYP)

c11 = [minimum(XYP[:,1]), maximum(XYP[:,2])]
cmn = [maximum(XYP[:,1]), minimum(XYP[:,2])]

h = Array{Plots.Plot{Plots.GRBackend}}(undef, J)
for i in 1:J
  h[i] = heatmap(LinRange(-2.0, 2.0, size(G,1)), LinRange(-2.0, 2.0, size(G, 2)), 
                G[:,:,i], c = :seismic, colorbar = false)
  plot!(h[i], XY[:,1], XY[:,2], c = :black, linestyle = :dash)
  plot!(h[i], c11[[1, 2, 2, 1, 1]], cmn[[1, 1, 2, 2, 1]], c = :black)
end

plot(h[1:5]...)
png("fried_eggs")

