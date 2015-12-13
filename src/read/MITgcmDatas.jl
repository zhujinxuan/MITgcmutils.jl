include("abstract.jl")


type MITgcmDatas <: DataFiles
  Dir :: ASCIIString 
  nx :: Int64
  ny :: Int64
  nz :: Int64
  dxspacing :: Float64
  dyspacing :: Float64
  dt :: Float64
  files :: Array{ASCIIString,1}
end


function MITgcmDatas(Dir; nx=200,ny=200, nz=60,
  dxspacing = 100.0, dyspacing = 100.0, dt = 30.0)
  files = map(ASCIIString,readdir(Dir))
  return MITgcmDatas(Dir, nx,ny,nz, dxspacing, dyspacing, dt,files)
end


export MITgcmDatas
export bread


