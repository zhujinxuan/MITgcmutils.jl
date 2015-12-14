abstract DimensionKmethod

type OnePointLR <: DimensionKmethod 
end
type OnePointR <: DimensionKmethod
end
type OnePointL <: DimensionKmethod
end
export DimensionKmethod, OnePointR, OnePointLR, OnePointL

function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointLR)

  pidx = map( x->:, size(u))
  df = zeros(size(u))

  midx = copy(pidx); midx[varidx] = 2:(size(u, varidx)-1)
  Lidx = copy(pidx); Lidx[varidx] = 1:(size(u, varidx)-2)
  Ridx = copy(pidx); Ridx[varidx] = 3:size(u, varidx)
  df[midx...] = (u[Ridx...] - u[Lidx...])/2

  midx = copy(pidx); midx[varidx] = 1
  Lidx = copy(pidx); Lidx[varidx] = 1
  Ridx = copy(pidx); Ridx[varidx] = 2
  df[midx...] = u[Ridx...] - u[Lidx...]

  midx = copy(pidx); midx[varidx] = size(u, varidx)
  Lidx = copy(pidx); Lidx[varidx] = size(u, varidx) -1
  Ridx = copy(pidx); Ridx[varidx] = size(u, varidx) 
  df[midx...] = u[Ridx...] - u[Lidx...]

  return df
end


function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointL)

  pidx = map( x->:, size(u))
  df = zeros(size(u))

  midx = copy(pidx); midx[varidx] = 2:(size(u, varidx))
  Lidx = copy(pidx); Lidx[varidx] = 1:(size(u, varidx)-1)
  Ridx = copy(pidx); Ridx[varidx] = 2:(size(u, varidx))
  df[midx...] = (u[Ridx...] - u[Lidx...])

  midx = copy(pidx); midx[varidx] = 1
  Lidx = copy(pidx); Lidx[varidx] = 1
  Ridx = copy(pidx); Ridx[varidx] = 2
  df[midx...] = u[Ridx...] - u[Lidx...]
  return df
end


function dimDiff( u :: Array{Float64}, varidx :: Int64, m :: OnePointR)

  pidx = map( x->:, size(u))
  df = zeros(size(u))

  midx = copy(pidx); midx[varidx] = 1:(size(u, varidx)-1)
  Lidx = copy(pidx); Lidx[varidx] = 1:(size(u, varidx)-1)
  Ridx = copy(pidx); Ridx[varidx] = 2:(size(u, varidx))
  df[midx...] = (u[Ridx...] - u[Lidx...])

  midx = copy(pidx); midx[varidx] = size(u, varidx)
  Lidx = copy(pidx); Lidx[varidx] = size(u, varidx) -1
  Ridx = copy(pidx); Ridx[varidx] = size(u, varidx) 
  df[midx...] = u[Ridx...] - u[Lidx...]
  return df
end


