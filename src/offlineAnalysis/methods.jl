function pv( u :: Array{Float64,2}, v :: Array{Float64,2}, m :: MITgcmDatas)
  dv = zeros(size(u))
  du = zeros(size(u))
  dxspacing = m.dxspacing
  dyspacing = m.dyspacing
  dv[2:end-1,:] = (v[3:end,:] - v[1:end-2,:])/(2*dxspacing)
  dv[1,:] = (v[2,:] - v[1,:])/dxspacing
  dv[end,:] = (v[end,:] - v[end-1,:])/dxspacing

  du[:,2:end-1] = (u[:,3:end] - u[:,1:end-2])/(2*dyspacing)
  du[:,1] = (u[:,2] - u[:,1])/dyspacing
  du[end,:] = (u[:,end] - u[:,end-1])/dyspacing
  return du+dv
end

function pv( u :: Array{Float64,3}, v :: Array{Float64,3}, m :: MITgcmDatas)
  res = zeros(size(u))
  for ii = 1:size(res,3)
    res[:,:,ii] = pv(u[:,:,ii],v[:,:,ii],m)
  end
  return res
end
export pv

