function verticalA( w :: Array{Float64,3}, T :: Array{Float64,3}, wintercept :: OnePointR = OnePointR())
  w1 = dimAverage(w,3, wintercept)
  return T.* w1
end


export verticalA
