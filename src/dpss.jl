
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

""" Eigenvalues/concentrations for the  Slepian sequences.   
They are computed as given in Percival 390, BUT THERE IS A TYPO IN PERCIVAL 390 """
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

