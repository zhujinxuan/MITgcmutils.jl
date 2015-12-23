function days(times :: Array{Int64,1}, m :: MITgcmDatas 
              ;min = 15, max = 20 , dt = m.dt )
  return (min*86400 .< (times * dt).< max*86400)
end

export days

function filestimes(m:: MITgcmDatas, varfile :: ASCIIString)
  rr = Regex("^$(varfile)\\.0.*\\.data\$")
  filenames = filter( rr, m.files)
  nsteps  = map(filenames) do x
    rr1 = Regex("^$(varfile)\\.0")
    x1 = replace(x, rr1, s"")
    p = match(r"[0-9]+",x1)
    parse(Int64,p.match)
  end
  return (filenames, nsteps )
end

function exactTimeStep(times :: Array{Int64,1}, m :: MITgcmDatas;
                       minstep :: Int64 = 56880, maxstep :: Int64 = 57000)
  return  minstep .<= times .<= maxstep
end
export exactTimeStep
