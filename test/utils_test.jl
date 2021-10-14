
using Random
Random.seed!(123)

N = 1024
K = 7
NW = 4.0
W = 0.08

ev, v = dpss_tapers(N, W, K, :both)

@testset "dpss test" begin
  
  # Eigenvalues
  @test ev[1:3] ≈ [0.031031839127290427, 0.03103311455746982, 0.031034387523034522]

  # Eigenvectors
  @test v[1:3] ≈[0.1588842745838488, 0.0011144647659785548, 1.2600650734933855e-6]

end


t = sort(randn(1024))
ev, v = gpss_orth(W, K, t, 0.1)

@testset "gpss test" begin 
  
  # Eigenvalues
  @test ev[1:3] ≈ [0.8488638017599491, 0.4569294729712821, 0.2890809140379784]

  # Eigenvectors
  @test v[1:3] ≈ [0.14268716855395974, 0.09159710145322413, 0.1682672984982918]

end

sleps = customsleps(5, 10.0, (16, 8), exact = true)

sleps_hodlr = customsleps(5, 10.0, (16, 8), exact = false)

@testset "2D Spherical Slepians" begin

  @test sleps[1][1:3] ≈ Complex{Float64}[1.4165175331754862 + 0.0im, 1.302170654345073 + 0.0im, 1.1413498706076248 + 0.0im]

  # Exact 2D Slepians
  @test sleps[2][2][1:3] ≈ Complex{Float64}[-0.10725724727346368 + 0.0im, -0.17892535707531507 + 0.0im, -0.23728562389432664 + 0.0im]

  # HODLR version
  @test sleps_hodlr[1][1:3] ≈ Complex{Float64}[1.4165175332530344 + 0.0im, 1.3021706541792697 + 0.0im, 1.1413498706671528 + 0.0im]

  @test sleps_hodlr[2][2][1:3] ≈ Complex{Float64}[-0.10725724727179216 + 0.0im, -0.1789253508692578 + 0.0im, -0.23728562378174453 + 0.0im]

end
