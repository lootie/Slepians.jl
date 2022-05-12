# Run the custom slepians - small modifications made to the main code:
# https://github.com/lootie/Slepians.jl/blob/master/src/dDimSleps.jl

function customsleps_ext(M::Int64, Kp, szs; prec = 1.0e-8, exact = false, lvl = 6, maxrank = 0, 
    no = nothing, sqwt = nothing, int = nothing, ev = nothing)
    # Create the kernel matrix, exactly or via HODLR approximation:
    _K    = KernelMatrix(no, no, Kp, dfn)
    K     = exact ? full(_K) : RHMatrix.rhodlr(_K, lvl, prec, maxrank)

    # Solve the eigenvalue problem (86) to obtain the slepians at the quad nodes:
     s     = eigsolve(z -> sqwt .* (K * (sqwt .* z)), length(sqwt), M, :LM, 
                        issymmetric = true)
    outputsize = (int == nothing) ? szs : int
    # Prepare even grid points, compute Slepians at those points:
    evpts = (ev != nothing) ? ev : vec(collect(product([range(-1.0, 1.0, length = s) for s in outputsize]...)))

    _K2   = KernelMatrix(evpts, no, Kp, dfn)
    K2    = exact ? full(_K2) : RHMatrix.rhodlr(_K2, lvl, prec, maxrank)
    sleps = [(K2 * (sqwt .* x[2])) ./ x[1] for x in zip(s[1], s[2])]
    return s[1], [reshape(sp, szs...) for sp in sleps]; #szs or outputsize?
end