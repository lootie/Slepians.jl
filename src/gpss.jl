
"""
    gpss(w, k, t, f; <keyword arguments>)

Generalized prolate spheroidal sequences on an unequal grid

...
# Arguments

## Positional Arguments

 - `w::Float64`: the bandwidth

 - `k::Int64`: number of Slepian tapers, must be <=2*bw*length(x) 

 - `t::Vector{Int64}`: vector containing the time indices

 - `f::Float64`: frequency at which the tapers are to be computed

## Keyword Arguments

 - `beta::Float64 = 0.5`: analysis half-bandwidth (similar to Nyquist rate)

...

...

# Outputs

 - `lambda::Vector{Float64}` the concentrations of the generalized prolate spheroidal
sequences

 - `u::Matrix{Float64}` the matrix containing the sequences themselves

 - `R` the Cholesky factor for the generalized eigenvalue problem

...

See also: [`gpss_orth`](@ref)

"""
function gpss(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}}, 
        f::Float64; beta::Float64 = 0.5)
  n           = length(t)
  a           = 2*w*ones(n,n)
  b = 2.0*beta*ones(n,n) .+ 0.0im
  for i = 1:n
      for j in (i+1):n
          a[i,j]  = sin.(2*pi*w*(t[i] .- t[j]))./(pi*(t[i] .- t[j]))
          b[i,j]  = exp.(-2*pi*1.0im*f*(t[i] .- t[j])).*
                      sin.(2*pi*beta*(t[i] .- t[j]))./(pi*(t[i] .- t[j])) 
      end
  end
  # If the Cholesky factorization fails, add a small number to the diagonal and
  # try again. Then remove that same value from the eigenvalues. 
  R, fudge = try (cholesky(Hermitian(b)), 0.0)
  catch
    (cholesky(Hermitian(b) + Matrix(I, size(b)...)*(1e-10)), 1.0)
  end
  V = eigen(Symmetric(real.(inv(R.L)*Symmetric(a)*inv(R.U))))
  lambda = V.values[end:-1:(end-k+1)] .- fudge*(1e-10)
  v = inv(R.U)*V.vectors[:,end:-1:(end-k+1)]
  u = copy(v)
  for i = 1:2:k
      if mean(real.(u[:,i])) < 0 
        u[:,i] = -u[:,i] 
      end
  end
  for i = 2:2:k-1
      if real(u[2,i] - u[1,i]) < 0
        u[:,i] = -u[:,i] 
      end
  end
  return (lambda, u, R)
end

"""
    gpss_orth(w, k, t, f; <keyword arguments>)

Generalized, orthogonalized prolate spheroidal sequences on an unequal grid

...
# Arguments

## Positional Arguments

 - `w::Float64`: the bandwidth

 - `k::Int64`: number of Slepian tapers, must be <=2*bw*length(x) 

 - `t::Vector{Int64}`: vector containing the time indices

 - `f::Float64`: frequency at which the tapers are to be computed

## Keyword Arguments

 - `beta::Float64 = 0.5`: analysis half-bandwidth (similar to Nyquist rate)

...

...

# Outputs

 - `lambda::Vector{Float64}` the concentrations of the generalized prolate spheroidal
sequences

 - `u::Matrix{Float64}` the matrix containing the sequences themselves, equivalent to 
u*R for the ordinary `gpss` routine.

...

See also: [`gpss`](@ref)

"""
function gpss_orth(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}}, 
        f::Float64; beta::Float64 = 0.5)
  n           = length(t)
  a           = 2*w*ones(n,n)
  b = 2.0*beta*ones(n,n) .+ 0.0im
  for i = 1:n
      for j in (i+1):n
          a[i,j]  = sin.(2*pi*w*(t[i] .- t[j]))./(pi*(t[i] .- t[j]))
          b[i,j]  = exp.(-2*pi*1.0im*f*(t[i] .- t[j])).*
                      sin.(2*pi*beta*(t[i] .- t[j]))./(pi*(t[i] .- t[j])) 
      end
  end
  # If the Cholesky factorization fails, add a small number to the diagonal and
  # try again. Then remove that same value from the eigenvalues. 
  R, fudge = try (cholesky(Hermitian(b)), 0.0)
  catch
    (cholesky(Hermitian(b) + Matrix(I, size(b)...)*(1e-10)), 1.0)
  end
  V = eigen(Symmetric(real.(inv(R.L)*Symmetric(a)*inv(R.U))))
  lambda = V.values[end:-1:(end-k+1)] .- fudge*(1e-10)
  u = V.vectors[:,end:-1:(end-k+1)]
  for i = 1:2:k
      if mean(real.(u[:,i])) < 0 
        u[:,i] = -u[:,i] 
      end
  end
  for i = 2:2:k-1
      if real(u[2,i] - u[1,i]) < 0
        u[:,i] = -u[:,i] 
      end
  end
  return (lambda, u)
end

