
N = 1024
K = 7
NW = 4.0

@testset "dpss test" begin 
  
  # Expected jackknife variance
  @test sin(pi/2) ≈ 1.0 

end
