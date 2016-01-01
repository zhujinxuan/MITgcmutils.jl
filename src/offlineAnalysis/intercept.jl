
abstract Position

type GridPosition
  x :: Float64
  y :: Float64
end
uGrid = GridPosition(0.0,0.5)
vGrid = GridPosition(0.5,0.0)
mGrid = GridPosition(0.0,0.0)
sGrid = GridPosition(0.5,0.5)


type GridMoving
  xshift :: Float64
  yshift :: Float64
  xband :: Float64
  yband :: Float64
end

function GridMoving( o :: GridPosition, t :: GridPosition, xband :: Float64, yband :: Float64)
  xshift = t.x - o.x 
  yshift = t.y - o.y 
  return GridMoving(xshift, yshift, xband, yband)
end


function intercept( ss :: Array{Float64,2}, g :: GridMoving; keep_when_notmoving :: Bool = false)
  signal = ss
  xband = g.xband
  yband = g.yband
  xshift = g.xshift
  yshift = g.yshift
  kx = collect(1:(size(signal,1)/xband)) / size(signal,1)
  ky = collect(1:(size(signal,2)/yband)) / size(signal,2)

  fresx = zeros(signal)
  if ((xshift == 0.0) & keep_when_notmoving)
    fresx = signal
  else
    for ix = kx
      psix = ix*pi*collect(1:size(signal,1))
      psix1 = ix*pi*(xshift + collect(1:size(signal,1)))
      fsinx  = mean(sin(psix) .* signal,1).* sin(psix1)/(0.5)
      fcosx  = mean(cos(psix) .* signal,1).* cos(psix1)/(0.5)
      fresx += fsinx + fcosx
    end
  end

  signal = fresx
  fresy = zeros(size(signal))
  if ((yshift == 0.0) & keep_when_notmoving)
    fresy = signal
  else
    for iy = ky
      psiy = iy*pi*collect(1:size(signal,2))
      psiy1 = iy*pi*(yshift + collect(1:size(signal,2)))
      psiy = psiy'
      psiy1 = psiy1'
      fsiny  = mean(sin(psiy) .* signal,2).* sin(psiy1)/(0.5)
      fcosy  = mean(cos(psiy) .* signal,2).* cos(psiy1)/(0.5)
      fresy += fsiny + fcosy
    end
  end

  return fresy
end

function intercept(ss:: Array{Float64,3}, g:: GridMoving; kws...)
  res = zeros(ss)
  for iz = 1:size(ss,3)
    res[:,:,iz] = intercept(ss[:,:,iz],g;kws...)
  end
  return res
end
export intercept
export GridMoving
