var documenterSearchIndex = {"docs":
[{"location":"#Slepians.jl","page":"Slepians.jl","title":"Slepians.jl","text":"","category":"section"},{"location":"","page":"Slepians.jl","title":"Slepians.jl","text":"This repo includes codes to compute discrete prolate spheroidal sequences, generalized prolate spheroidal sequences (unequal sampling), and solves the special case of d-dimensional concentration in the Cartesian plane where the function is restricted to the disc in the spectral domain.","category":"page"},{"location":"","page":"Slepians.jl","title":"Slepians.jl","text":"These docs are under construction, but I will include the docstrings for the main functions here.","category":"page"},{"location":"#Index-of-functions","page":"Slepians.jl","title":"Index of functions","text":"","category":"section"},{"location":"","page":"Slepians.jl","title":"Slepians.jl","text":"Modules = [Slepians]\nOrder   = [:function, :type]","category":"page"},{"location":"#Slepians.blob","page":"Slepians.jl","title":"Slepians.blob","text":"blob(N, Nj)\n\nMakes moving picture of a random blob by superposition of random circles\n\nArguments\n\nN number of loops for, and if movie, default = 100\nNj smoothness, roughly (?)\n\nOutputs\n\nx\ny \n\n\n\n\n\n","category":"function"},{"location":"#Slepians.closedcurve_2d-Tuple{Any, Any}","page":"Slepians.jl","title":"Slepians.closedcurve_2d","text":"closedcurve_2d(xs, ys)\n\nGenerate a closed curve based on (x, y) points given and Euclidean distance\n\nArguments\n\nxs x coordinates of the points\nys y coordinates of the points\n\nOutputs\n\nxpath x coordinates of the path\nypath y coordinates of the path\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.conv-Tuple{Any, Any}","page":"Slepians.jl","title":"Slepians.conv","text":"conv(x,y)\n\nConvolves two arrays of the same length using FFT algorithm\n\nArguments\n\nx::Array{Number}: first array\ny::Array{Number}: second array\n\nOutputs\n\nArray containing the real part of the convolution of x and y \n\n\n\n\n\n","category":"method"},{"location":"#Slepians.degamini-Tuple{Any}","page":"Slepians.jl","title":"Slepians.degamini","text":"degamini(v)\n\nArguments\n\nv A vector with repeated entries\n\nOutputs\n\ndv The same vector with all the repeat sset to 1\nfoldi A same-dimensional vector with how many repeats there are\nbe A matrix with begin and end indices into the original vector\n\nExample\n\ndv, foldi, be = degamini([1, 2, 2, 3])\n\nanswer is dv = [1, 2, 3]; foldi = [1, 2, 1]; be = [1 1; 2 3; 4 4]\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.dpss_eigval-NTuple{4, Any}","page":"Slepians.jl","title":"Slepians.dpss_eigval","text":"dpss_eigval(dpVecs, n, nw, ntapers)\n\nEigenvalues/concentrations for the Slepian sequences, given the vectors (Percival 390)\n\nArguments\n\ndpVecs: (n x ntapers) Matrix in which Slepian eigenvectors are the columns\nn: Integer length of the data\nnw: Float time bandwidth product\nntapers: Integer number of tapers\n\nOutputs\n\neeigvalss: Vector of lengthh ntapers containing the eigenvalues\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.dpss_tapers","page":"Slepians.jl","title":"Slepians.dpss_tapers","text":"dpss_tapers(n,w,k,tap_or_egval)\n\nSimply compute discrete prolate spheroidal sequence tapers, eigenvalues \n\n...\n\nArguments\n\nn::Int64: Length of the taper\nnw::Float64: Time-bandwidth product\nk::Int64: Number of tapers\ntap_or_egval::Symbol = :tap: Either :tap, :egval, or :both\n\n...\n\n...\n\nOutputs\n\nvv::Vector{Float64}: The matrix of eigenvalues, if taporegval is set to :tap\ndpss_eigval: Struct conaining the dpss tapers\n\n...\n\n\n\n\n\n","category":"function"},{"location":"#Slepians.findnst-Tuple{Vector{T} where T, Vector{T} where T, Any, Any}","page":"Slepians.jl","title":"Slepians.findnst","text":"findnst(xs, ys, x0, y0)\n\nFind the nearest 2D point in arrays xs, ys to (x0, y0) Euclidean distance \n\n\n\n\n\n","category":"method"},{"location":"#Slepians.gamini","page":"Slepians.jl","title":"Slepians.gamini","text":"gamini(data, folding)\n\nArguments\n\n-data: some data vector -folding: The replication factor for every element of the data; if a scalar thsi applies to all of the elements (default 3), if zero or negative, no replication occurs\n\nOutputs\n\n-bigger\n\nExample usage\n\na, b = degamini(gamini([1, 2, 3, 1, 4, 5], [1, 2, 3, 2, 4, 2]))\n\nOne gets [1, 2, 2, 3, 3, 3, 1, 1, 4, 4, 4, 4, 5, 5] as the intermediate result.\n\nSee also\n\n@degamini, @gamini2\n\n\n\n\n\n","category":"function"},{"location":"#Slepians.get_plan-Tuple{Any}","page":"Slepians.jl","title":"Slepians.get_plan","text":"get_plan(n)\n\nObtain an FFT plan\n\nArguments\n\nn::Int64: Length of the FFT plan\n\nOutputs\n\nTuple containing two ComplexF64 arrays of length n and an FFT plan\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.getbdypts_2d-Union{Tuple{Matrix{T}}, Tuple{T}} where T<:Number","page":"Slepians.jl","title":"Slepians.getbdypts_2d","text":"getbdypts_2d(mask)\n\nStarting with a matrix of zeros and ones, get the points on the boundary\n\nArguments\n\n- `mask` a matrix of zeros and ones\n\nOutputs\n\n- `x` x-coordinates of the boundary points\n- `y` y-coordinates of the boundary points\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.gpss-Tuple{Float64, Int64, Union{Vector{Float64}, Vector{Int64}}, Float64}","page":"Slepians.jl","title":"Slepians.gpss","text":"gpss(w, k, t, f; <keyword arguments>)\n\nGeneralized prolate spheroidal sequences on an unequal grid\n\n...\n\nArguments\n\nPositional Arguments\n\nw::Float64: the bandwidth\nk::Int64: number of Slepian tapers, must be <=2bwlength(x) \nt::Vector{Int64}: vector containing the time indices\nf::Float64: frequency at which the tapers are to be computed\n\nKeyword Arguments\n\nbeta::Float64 = 0.5: analysis half-bandwidth (similar to Nyquist rate)\n\n...\n\n...\n\nOutputs\n\nlambda::Vector{Float64} the concentrations of the generalized prolate spheroidal\n\nsequences\n\nu::Matrix{Float64} the matrix containing the sequences themselves\nR the Cholesky factor for the generalized eigenvalue problem\n\n...\n\nSee also: gpss_orth\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.gpss_orth-Tuple{Float64, Int64, Union{Vector{Float64}, Vector{Int64}}, Float64}","page":"Slepians.jl","title":"Slepians.gpss_orth","text":"gpss_orth(w, k, t, f; <keyword arguments>)\n\nGeneralized, orthogonalized prolate spheroidal sequences on an unequal grid\n\n...\n\nArguments\n\nPositional Arguments\n\nw::Float64: the bandwidth\nk::Int64: number of Slepian tapers, must be <=2bwlength(x) \nt::Vector{Int64}: vector containing the time indices\nf::Float64: frequency at which the tapers are to be computed\n\nKeyword Arguments\n\nbeta::Float64 = 0.5: analysis half-bandwidth (similar to Nyquist rate)\n\n...\n\n...\n\nOutputs\n\nlambda::Vector{Float64} the concentrations of the generalized prolate spheroidal\n\nsequences\n\nu::Matrix{Float64} the matrix containing the sequences themselves, equivalent to \n\nu*R for the ordinary gpss routine.\n\n...\n\nSee also: gpss\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.interp1","page":"Slepians.jl","title":"Slepians.interp1","text":"interp1()\n\n1D data interpolation, linear https://www.mathworks.com/help/matlab/ref/interp1.html?stid=docta\n\nArguments\n\nVector x contains the sample points, and \nv contains the corresponding values, v(x). \nVector xq contains the coordinates of the query points.\n\nOutputs\n\ninterpolated values of a 1-D function at specific query points using linear interpolation\n\nExample usage\n\ninterp1(LinRange(0, 1, 10), log10.(LinRange(0, 1, 10)), [0.5, 0.9], :linear)\n\n\n\n\n\n","category":"function"},{"location":"#Slepians.interp2","page":"Slepians.jl","title":"Slepians.interp2","text":"interp2(z, xy, z0)\n\nArguments\n\nz::Vector vector of z's\nyx::Matrix y, x coordinates at each of the z's\nz0<:Number level at which to interpolate z0\n\nOutputs\n\nzyx coordinates corresponding to level z0\n\n\n\n\n\n","category":"function"},{"location":"#Slepians.interpcontour-NTuple{4, Any}","page":"Slepians.jl","title":"Slepians.interpcontour","text":"interpcontour(z, z0, thph, N)\n\nArguments\n\nz::Vector the two z-values between which the level z0 falls. \nz0<:Number the z-level of the desired contour\nthph::Matrix the zyx coordinates (first N are at level z[1] second are at level z[2])\nN::Int64 number of points for each contour\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.matranges-Tuple{Any}","page":"Slepians.jl","title":"Slepians.matranges","text":"matranges(ranges)\n\nMakes an index vector with monotonically increasing indices between pairs of numbers supplied as input. From slepian_alpha\n\nArguments\n\nranges\nOutputs\n\nExample\n\nmatranges([1, 4, 1, 2, -1, 2])\n# answer is [1, 2, 3, 4, 1, 2, -1, 0, 1, 2]\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.mdslepian-Tuple{Any, Any, Any}","page":"Slepians.jl","title":"Slepians.mdslepian","text":"mdslepian(w, k, t)\n\nGeneralized prolate spheroidal sequences for the 1D missing data problem\n\n...\n\nArguments\n\nPositional Arguments\n\nw::Float64: the bandwidth\nk::Int64: number of Slepian tapers, must be <=2bwlength(x) \nt::Vector{Int64}: vector containing the time indices\n\n...\n\n...\n\nOutputs\n\nlambda,u::Tuple{Vector{Float64}, Vector{Float64}}: tuple containing the \n\nconcentrations and the tapers\n\n...\n\nSee also: mdmultispec, gpss\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.phicurve-Tuple{Any, Any}","page":"Slepians.jl","title":"Slepians.phicurve","text":"phicurve(thph, th)\n\nAdapted from slepian_alpha; finds the longitude crossings and thus integration domains of a closed curve parameterized in colatitude/longitude space at certian query points of colatitude.\n\nArguments\n\nthph:: Colatitude/longitude of the closed curve (degrees)\nth:: Colatitude at which crossings are required (degrees)\n\nOutputs\n\nphint A matrix with crossings/intervals and zeros of dimensions MxN where M = length(th) and N can be anything depending on the oscillations of the curve\nthp Colatitude matrix for hatched plotting, if possible\nphp Longitude matrix for hatched plotting, if possible\nforreal Indices of the ones that are real (could be at zero)\n\nDepends on: @sub2ind, @interp1, @degamini, @matranges, and possibly @blob (demo2)\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.quadpts-NTuple{7, Any}","page":"Slepians.jl","title":"Slepians.quadpts","text":"quadpts(qx, Nqx, qy, forreal, xints, wqx, wqy)\n\nScale the quaduature points and weights \n\nArguments\n\nqx \nNqx\nqy\nforreal\nxints\nwqx\nwqy\n\nOutputs\n\nQX Quadrature points in x\nQY Quadrature points in y\nw Quadrature weights\nNrun Number of horizontal line segments in the domain\n\n\n\n\n\n","category":"method"},{"location":"#Slepians.randcirc","page":"Slepians.jl","title":"Slepians.randcirc","text":"randcirc(xm, ym, r,dr, N)\n\nArguments\n\nxm horizontal positon of the center\nym vertical position of the center\nr radius\ndr size of random perturbations around the radius\nN number of random spike points\n\nOutputs\n\nx x-coordinate\ny y-coordinate\n\nRelated\n\n@blob\n\n\n\n\n\n","category":"function"},{"location":"#Slepians.sub2ind-Tuple{Any, Any, Any}","page":"Slepians.jl","title":"Slepians.sub2ind","text":"sub2ind(A, row, col)\n\nConvert to linear indices  https://www.mathworks.com/help/matlab/ref/sub2ind.html\n\n\n\n\n\n","category":"method"}]
}