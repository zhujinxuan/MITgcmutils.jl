function days(times :: Array{Float64,1};min = 15, max = 20 )
  return (min*86400 .< times .< max*86400)
end

export days

function filestimes(m:: MITgcmDatas, varfile :: ASCIIString)
  filenames = filter( Regex("^$(varfile).*\\.data\$"), m.files)
  nsteps  = map(filenames) do x
    p = match(r"[0-9]+",x)
    parse(Int64,p.match)
  end
  return (filenames, nsteps * m.dt)
end


