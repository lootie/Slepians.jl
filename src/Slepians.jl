module Slepians

  using FastGaussQuadrature, IterTools, Statistics 
  using FFTW, LinearAlgebra, SpecialFunctions, Arpack, KrylovKit, KernelMatrices 
  using Interpolations, PCHIPInterpolation


  include("RHMatrix.jl")
  using .RHMatrix

  include("dpss.jl")
  export dpss_tapers

  include("gpss.jl")
  export gpss, gpss_orth, mdslepian

  include("dDimSleps.jl")
  export dfn, efn, givewts, getnodeswts, customsleps 

  include("relevant_slepian_alpha.jl")
  export interp1, gamini, degamini, matranges, randcirc, blob, phicurve, quadpts
  export get_quadrature_nodes_2D

  include("bdy2d.jl")
  export getbdypts_2d, findnst, closedcurve_2d, mask2closedcurve

  include("ext.jl")
  export customsleps_ext

end 

# module
