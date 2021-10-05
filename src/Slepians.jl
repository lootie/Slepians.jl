module Slepians

  using FastGaussQuadrature, IterTools 
  using FFTW, LinearAlgebra, SpecialFuncions, Arpack, KrylovKit, KernelMatrices 

  include("RHmatrix.jl")

  include("dpss.jl")
  export dpss_tapers

  include("gpss.jl")
  export gpss, gpss_orth

end # module
