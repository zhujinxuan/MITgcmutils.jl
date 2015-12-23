function vor( u :: Array{Float64,2}, v :: Array{Float64,2}, m :: MITgcmDatas, partial :: OnePointR = OnePointR())
  dudy = dimDiff(u, 2, partial) /m.dyspacing
  dvdx = dimDiff(v, 1, partial) /m.dxspacing
  return (dvdx - dudy)
end

function vor( u :: Array{Float64,2}, v :: Array{Float64,2}, m :: MITgcmDatas, partial :: DimensionKmethod )
  dudy = dimDiff(u, 2, partial) /m.dyspacing
  dvdx = dimDiff(v, 1, partial) /m.dxspacing
  return (dvdx - dudy)
end

function MITgcmVor( u :: Array{Float64,3}, v :: Array{Float64,3}, m :: MITgcmDatas; dzspacing :: Float64 = 2.5, partial :: OnePointR = OnePointR())

  dudy = dimDiff(u, 2, partial) /m.dyspacing
  dvdx = dimDiff(v, 1, partial) /m.dxspacing
  dudz = - dimDiff(u, 3, partial) /  dzspacing
  dvdz = - dimDiff(v, 3, partial) /  dzspacing
  return ( - dvdz, dudz, (dvdx - dudy) )
end

function MITgcmPV( u :: Array{Float64,3}, v :: Array{Float64,3},
             T :: Array{Float64,3}, m :: MITgcmDatas; 
             dzspacing :: Float64 = 2.5, f :: Float64 = 7.29e-5)

  vorA  = begin  
    vv = MITgcmVor(u,v,m, partial = OnePointR())
    (vv[1], vv[2], vv[3]+f)
  end

  dT = begin
    TT = map( x->dimDiff(T,x, OnePointLR()), (1,2,3))
    (TT[1]/m.dxspacing, TT[2]/m.dyspacing, - TT[3]/dzspacing)
  end

  pvres = mapreduce(+, zip(vorA, dT) ) do x
    x[1] .* x[2]
  end
  return pvres
end

export vor
export MITgcmPV, MITgcmVor

