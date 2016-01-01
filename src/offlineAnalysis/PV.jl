function vor( u :: Array{Float64,2}, v :: Array{Float64,2}, m :: MITgcmDatas)
  dudy = dimDiff(u, 2, OnePointR()) /m.dyspacing
  dvdx = dimDiff(v, 1, OnePointR()) /m.dxspacing
  return (dvdx - dudy)
end

function vor( u :: Array{Float64,3}, v :: Array{Float64,3}, m :: MITgcmDatas)
  vv = zeros(size(u))
  for iz = 1: size(u,3)
    vv[:,:,iz] = vor(u[:,:,iz],v[:,:,iz],m)
  end
  return vv 
end


#= include("intercept.jl") =#
function MITgcmPV( u :: Array{Float64,3}, v :: Array{Float64,3},
             T :: Array{Float64,3}, m :: MITgcmDatas; 
             dzspacing :: Float64 = 2.5, f :: Float64 = 7.29e-5, BbyT :: Float64 = 2.0e-4 * 9.8)
  parIndent = OnePointR()
  dudy = begin
    fdy =  dimDiff(u, 2, OnePointLR()) /(m.dyspacing)
    dimAverage(fdy,1, parIndent)
    dimMoves(fdy,2,OnePointR())
  end
  dvdx = begin
    fdx = dimDiff(v,1,OnePointLR())/(m.dxspacing)
    dimAverage(fdx, 2,parIndent)
    dimMoves(fdx,1,OnePointR())
  end
  dudz = begin
    fdz = -dimDiff(u,3,OnePointLR())/dzspacing 
    dimAverage(fdz,1, parIndent)
    dimMoves(fdz,2,OnePointR())
  end
  dvdz = begin
    fdz = -dimDiff(v,3,OnePointLR())/dzspacing
    dimAverage(fdz,2, parIndent)
    dimMoves(fdz,1,OnePointR())
  end
  avor = ( -dvdz, dudz, f+(dvdx - dudy))


  dB = begin
    TT = map( x-> BbyT * dimDiff(T,x, OnePointLR()), (1,2,3))
    (TT[1]/m.dxspacing, TT[2]/m.dyspacing, -TT[3]/dzspacing)
  end

  return mapreduce(+, zip(avor, dB))do x
    x[1].*x[2]
  end
end

export vor
export MITgcmPV, MITgcmVor

include("intercept.jl")
function MFPV( u :: Array{Float64,3}, v :: Array{Float64,3},
             T :: Array{Float64,3}, m :: MITgcmDatas; 
             dzspacing :: Float64 = 2.5, f :: Float64 = 7.29e-5, BbyT :: Float64 = 2.0e-4 * 9.8,
             whether_smooth_T :: Bool = true
             )
  dudy = begin
    fdy =  dimDiff(u, 2, OnePointLR()) /(m.dyspacing)
    shift = GridMoving(uGrid, GridPosition(mGrid.x, mGrid.y - 0.5),2.0,2.0)
    intercept(fdy, shift)/(m.dyspacing)
  end
  dvdx = begin
    fdx = dimDiff(v,1,OnePointLR())/(m.dxspacing)
    shift = GridMoving(vGrid, GridPosition(mGrid.x - 0.5, mGrid.y),2.0,2.0) 
    intercept(fdx, shift)
  end
  dudz = begin
    fdz = -dimDiff(u,3,OnePointLR())/dzspacing 
    shift = GridMoving(uGrid, mGrid,2.,2.)
    intercept(fdz,shift)
  end
  dvdz = begin
    fdz = -dimDiff(v,3,OnePointLR())/dzspacing
    shift = GridMoving(vGrid, mGrid,2.,2.)
    intercept(fdz,shift)/dzspacing
  end
  avor = ( -dvdz, dudz, f+(dvdx - dudy))


  dB = begin
    T1 = intercept(T,GridMoving(0.0, 0.0, 2.0,2.0);  keep_when_notmoving=!whether_smooth_T)
    TT = map( x-> BbyT * dimDiff(T1,x, OnePointLR()), (1,2,3))
    (TT[1]/m.dxspacing, TT[2]/m.dyspacing, -TT[3]/dzspacing)
  end

  return mapreduce(+, zip(avor, dB))do x
    x[1].*x[2]
  end
end
export MFPV
