
using Random
Random.seed!(123)

N = 1024
K = 7
NW = 4.0
W = 0.08

t = sort(randn(1024))

ev, v = gpss_orth(W, K, t, 0.1)

@testset "gpss test" begin 
  
  # Expected jackknife variance
  @test ev[1] ≈ [0.8488638017599491, 0.4569294729712821, 0.2890809140379784]

  @test v[1:3] ≈ [0.13416225814078597, 0.09152484752991513, 0.1582872094843311]

end
