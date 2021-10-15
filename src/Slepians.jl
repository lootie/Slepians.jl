module Slepians

  using FastGaussQuadrature, IterTools, Statistics 
  using FFTW, LinearAlgebra, SpecialFunctions, Arpack, KrylovKit, KernelMatrices 

  include("RHMatrix.jl")
  using .RHMatrix

  include("dpss.jl")
  export dpss_tapers

  include("gpss.jl")
  export gpss, gpss_orth, mdslepian

  include("dDimSleps.jl")
  export dfn, efn, givewts, getnodeswts, customsleps 

end 

# module
