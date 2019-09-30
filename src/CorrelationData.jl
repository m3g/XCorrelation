
struct CorrelationData

  nCA :: Int64
  npairs :: Int64
  ncontacts :: Int64
  index :: Vector{Int64}
  dist :: Vector{Float64}
  bin :: Vector{Float64}
  resnum_from_seq :: Vector{Int64}
  ContactDistances :: Matrix{Float64}
  ContactBin :: Matrix{Float64}
  Contacts :: Matrix{Int64}
  ContactDref :: Vector{Float64}

end

