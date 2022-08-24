
# GNU GPL v2 licenced to C. Haley 10/2021
# This code originally belonged to Multitaper.jl, by Haley & Geoga

"""
    get_plan(n)

Obtain an FFT plan

# Arguments
- `n::Int64`: Length of the FFT plan

# Outputs
- Tuple containing two ComplexF64 arrays of length `n` and an FFT plan
"""
function get_plan(n)
  work1 = zeros(ComplexF64, n)
  work2 = zeros(ComplexF64, n)
  plan  = plan_fft!(work1)
  (work1, work2, plan)
end


"""
    conv(x,y)

Convolves two arrays of the same length using FFT algorithm

# Arguments
- `x::Array{Number}`: first array
- `y::Array{Number}`: second array

# Outputs
- Array containing the real part of the convolution of x and y 

"""
function conv(x,y)
  @assert length(x) == length(y) "This convolution only supports equal length"
  len = nextprod([2,3,5,7], 2*length(x)-1)
  (work1, work2, plan) = get_plan(len)
  fill!(work1, zero(ComplexF64))
  fill!(work2, zero(ComplexF64))
  work1[1:length(x)] .= x
  work2[1:length(y)] .= y
  plan*work1
  plan*work2
  work1.*=work2
  plan\work1
  isreal(x) ? real(work1) : work1
end

""" 

    dpss_eigval(dpVecs, n, nw, ntapers)

Eigenvalues/concentrations for the Slepian sequences, given the vectors (Percival 390)

# Arguments
- `dpVecs`: (n x ntapers) Matrix in which Slepian eigenvectors are the columns
- `n`: Integer length of the data
- `nw`: Float time bandwidth product
- `ntapers`: Integer number of tapers

# Outputs
- `eeigvalss`: Vector of lengthh `ntapers` containing the eigenvalues

"""
function dpss_eigval(dpVecs, n, nw, ntapers)
  sincvec  = [sin(2*pi*nw*j/n)/(pi*j) for j in 1:(n-1)]
  q_tau    = mapreduce(x -> conv(x, x)[n:(2*n-1)], hcat, eachcol(dpVecs))
  eigvalss = map(j->2*(nw/n + dot(q_tau[1,j].*q_tau[2:n,j], sincvec)), 1:ntapers)
  return eigvalss
end


""" 

    dpss_tapers(n,w,k,tap_or_egval)

Simply compute discrete prolate spheroidal sequence tapers, eigenvalues 

...

# Arguments

 - `n::Int64`: Length of the taper

 - `nw::Float64`: Time-bandwidth product

 - `k::Int64`: Number of tapers

 - `tap_or_egval::Symbol = :tap`: Either :tap, :egval, or :both
...

...

# Outputs

 - `vv::Vector{Float64}`: The matrix of eigenvalues, if tap_or_egval is set to :tap

 - `dpss_eigval`: Struct conaining the dpss tapers

...

"""
function dpss_tapers(n, nw, k, tap_or_egval::Symbol=:tap) 
  stdm  = SymTridiagonal([cos(2*pi*(nw/n))*abs2(0.5*(n-1)-(j-1)) for j in 1:n],
                          [0.5*j*(n-j) for j in 1:(n-1)])
  # Note that eigvals seems to sort ascending.
  vv = reverse(eigvecs(stdm, reverse(eigvals(stdm, (n-k+1):n))),dims=2)
  if tap_or_egval == :tap
    return vv  
  elseif tap_or_egval == :egval
    return dpss_eigval(vv, n, nw, k)
  elseif tap_or_egval == :both
    return(vv, dpss_eigval(vv, n, nw, k))
  end
end

