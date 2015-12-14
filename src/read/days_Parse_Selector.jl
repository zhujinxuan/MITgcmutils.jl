function days(times :: Array{Float64,1};min = 20, max = 30 )
  return (min*86400 .< times .< max*86400)
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
  return (filenames, nsteps * m.dt)
end


