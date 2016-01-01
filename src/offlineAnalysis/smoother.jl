abstract KernelFilter 
type AverageFilter <: KernelFilter
  xband :: Int64
  yband :: Int64
end

function AverageFilter()
  return AverageFilter(1,1)
end

type FourierFilter <: KernelFilter
  xband :: Float64
  yband :: Float64
end
function FourierFilter()
  return FourierFilter(2.0,2.0)
end

function xysmooth( signal :: Array{Float64,3}, filter :: KernelFilter)
  res = zeros(size(signal))
  for iz = 1:size(signal,3)
    res[:,:,iz] = xysmooth(signal[:,:,iz], filter)
  end
  return res
end

function xysmooth( x:: Array{Float64,2}, filter :: AverageFilter)
  y = zeros(size(x))
  xband = filter.xband
  yband = filter.yband
  for ix = 1:size(x,1)
    for iy = 1:size(x,2)
      px :: UnitRange{Int64} = max(1,ix-xband):min(size(x,1),ix + xband)
      py :: UnitRange{Int64} = max(1,iy-yband):min(size(x,2),iy + yband)
      y[ix,iy] = mean(x[px,py])
    end
  end
  return y
end


function xysmooth( ss :: Array{Float64,2}, filter :: FourierFilter)

  xband = filter.xband
  yband = filter.yband
  signal =ss

  kx = collect(1:(size(signal,1)/xband)) / size(signal,1)
  ky = collect(1:(size(signal,2)/yband)) / size(signal,2)


  fresx = zeros(size(signal))
  for ix = kx
    psix = ix*pi*collect(1:size(signal,1))
    fsinx  = mean(sin(psix) .* signal,1).* sin(psix)/(0.5)
    fcosx  = mean(cos(psix) .* signal,1).* cos(psix)/(0.5)
    fresx += fsinx + fcosx
  end
  signal = fresx

  fresy = zeros(size(signal))
  for iy = ky
    psiy = iy*pi*collect(1:size(signal,2))
    psiy = psiy'
    fsiny  = mean(sin(psiy) .* signal,2).* sin(psiy)/(0.5)
    fcosy  = mean(cos(psiy) .* signal,2).* cos(psiy)/(0.5)
    fresy += fsiny + fcosy
  end


  return fresy
end

export xysmooth
export AverageFilter, FourierFilter
export KernelFilter
