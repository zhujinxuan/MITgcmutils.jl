abstract KernelFilter 
type AverageFilter <: KernelFilter
  xband :: Int64
  yband :: Int64
end

function AverageFilter()
  return AverageFilter(1,1)
end


function xysmoothing( x:: Array{Float64,3}, filter :: AverageFilter)
  y = zeros(size(x))
  xband = filter.xband
  yband = filter.yband
  for ix = 1:size(x,1)
    for iy = 1:size(x,2)
      px = collect( max(1,ix-xband):min(size(x,1),ix + xband))
      py = collect( max(1,iy-yband):min(size(x,2),iy + yband))
      y[ix,iy,:] = mean(x[px,py,:],(1,2))
    end
  end
  return y
end

export xysmoothing
export AverageFilter
export KernelFilter
