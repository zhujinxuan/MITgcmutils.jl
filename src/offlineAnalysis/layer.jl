function layerMean( x :: Array{Float64, 3}, region = (1,2))
   return squeeze( mean(x, region), region)
end

export layerMean

function layermixing( x :: Array{Float64, 3}, cutoff :: Float64)
  h = zeros(Int64, size(x,1,2))
  for ix = 1:size(h,1)
    for iy = 1:size(h,2)
      h[ix,iy] = findlast( x[ix,iy,:]  .< cutoff)
    end
  end
  return h
end
export layermixing
