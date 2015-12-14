include("MITgcmDatas.jl")
include("days_Parse_Selector.jl")

using DataFrames
typealias Tsize Tuple{Int64, Vararg{Int64}}

function bread(singleFile :: ASCIIString, m :: MITgcmDatas; 
               varoffset :: Int64 = 1, 
               varsize :: Tsize = (m.nx, m.ny, m.nz) )
  open(singleFile) do fid
    seek(fid,8*(varoffset-1))
    x = read(fid, Float64, varsize)
    map(ntoh,x)
  end
end

function bread( m :: MITgcmDatas, varfile :: ASCIIString; 
                selector :: Function = days , 
                varoffset :: Int64 = 1, 
                whether_mean :: Bool = false, 
                varsize :: Tsize = (m.nx, m.ny, m.nz) )
  (files, times) = filestimes(m, varfile)
  (fs,ts) = map((files, times)) do x
    x[selector(times)]
  end
  println("Reading $(length(ts)) files")

  if (whether_mean)
    result = zeros(varsize...)
    for ii = 1:length(ts)
      result[:,:,:] +=  bread("$(m.Dir)/$(fs[ii])",m, varoffset = varoffset, varsize = varsize)
    end
    return result/length(ts)
  else
    result = zeros(varsize..., length(ts))
    for ii = 1:length(ts)
      result[:,:,:,ii] = bread("$(m.Dir)/$(fs[ii])",m, varoffset = varoffset, varsize = varsize)
    end
    return result
  end

end
export bread
