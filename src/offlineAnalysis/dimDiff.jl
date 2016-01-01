abstract DimensionKmethod

typealias TupleRanges Tuple{UnitRange{Int64},Vararg{UnitRange{Int64}}}

immutable OnePointLR <: DimensionKmethod end
immutable OnePointR <: DimensionKmethod end
immutable OnePointL <: DimensionKmethod end
export DimensionKmethod, OnePointR, OnePointLR, OnePointL

function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointLR)

  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx, Ridx) = map((2:(dimlength-1), 1:(dimlength -2), 3:(dimlength))) do  dimll
    pp=map(enumerate(size(u))) do xx
      (xx[1] == varidx) ? dimll : 1:xx[2]
    end
    (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])/2

  (midx, Lidx, Ridx) = map((1:1, 1:1, 2:2)) do  dimll
    pp =map(enumerate(size(u))) do xx
      (xx[1] == varidx) ? dimll : 1:xx[2]
    end
    (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  (midx, Lidx, Ridx) =
      map((dimlength:dimlength, (dimlength-1):(dimlength-1),dimlength:dimlength)) do dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  return df
end


function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointL)

  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx, Ridx) = 
      map((2:dimlength, 1:(dimlength -1), 2:dimlength)) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  (midx, Lidx, Ridx) = map((1:1, 1:1, 2:2)) do  dimll
    pp=map(enumerate(size(u))) do xx
      (xx[1] == varidx) ? dimll : 1:xx[2]
    end
    (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  return df
end


function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointR)
  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx, Ridx) = 
      map((1:(dimlength-1), 1:(dimlength -1), 2:(dimlength))) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : (1:xx[2])
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  (midx, Lidx, Ridx) =
      map((dimlength:dimlength, (dimlength-1):(dimlength-1),dimlength:dimlength)) do dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] - u[Lidx...])

  return df
end

export dimDiff
function dimAverage( u :: Array{Float64}, varidx :: Int64, m :: OnePointR)
  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx, Ridx) = 
      map((1:(dimlength-1), 1:(dimlength -1), 2:(dimlength))) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : (1:xx[2])
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] + u[Lidx...])/2.0

  (midx, Lidx, Ridx) =
      map((dimlength:dimlength, (dimlength-1):(dimlength-1),dimlength:dimlength)) do dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] + u[Lidx...])/2.0

  return df
end

function dimAverage( u :: Array{Float64}, varidx :: Int64, m :: OnePointL)
  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx, Ridx) = 
      map((2:dimlength, 1:(dimlength -1), 2:dimlength)) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : (1:xx[2])
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] + u[Lidx...])/2.0

  (midx, Lidx, Ridx) = map((1:1, 1:1, 2:2)) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = (u[Ridx...] + u[Lidx...])/2.0

  return df
end
export dimAverage

function dimMoves( u :: Array{Float64}, varidx :: Int64, m:: OnePointL)
  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Lidx) = 
      map((2:dimlength, 1:(dimlength -1))) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : (1:xx[2])
        end
        (pp...)
  end
  df[midx...] = u[Lidx...] 

  (midx, Lidx) = map((1:1, 1:1)) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = u[Lidx...] 

  return df
end

function dimMoves( u :: Array{Float64}, varidx :: Int64, m:: OnePointR)
  df = zeros(size(u))
  dimlength =size(u,varidx) 
  local midx :: TupleRanges
  local Lidx :: TupleRanges
  local Ridx :: TupleRanges

  (midx, Ridx) = 
      map((1:(dimlength-1), 2:(dimlength))) do  dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : (1:xx[2])
        end
        (pp...)
  end
  df[midx...] = u[Ridx...]

  (midx, Ridx) =
      map((dimlength:dimlength ,dimlength:dimlength)) do dimll
        pp = map(enumerate(size(u))) do xx
          (xx[1] == varidx) ? dimll : 1:xx[2]
        end
        (pp...)
  end
  df[midx...] = u[Ridx...] 
  return df
end
