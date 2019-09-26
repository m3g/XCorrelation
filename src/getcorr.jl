#
# Check if there is correlation data for these indices
#

function findi4(Correlation :: CorrelationData, i :: Int64, j :: Int64, k :: Int64, l :: Int64)
  i4 = index4D(Correlation.nCA,i,j,k,l)
  index = searchsortedfirst(Correlation.index,i4)
  if Correlation.index[index] == i4
    return index
  else
    return nothing
  end
end

#
# Convert indexes of PDB sequences in the indexes of sequential CAs
#

function iseq_from_resnum(Correlation :: CorrelationData, resnum :: Int64)
  iseq = findfirst( i -> i == resnum , Correlation.resnum_from_seq )
  if iseq == nothing
    error(" Could not find sequence number for resnum = ", resnum)
  end
  return iseq
end

#
# Get the corrleation of a specific pair of constraints defined by [i, j], [k, l]
#

function getcorr(Correlation :: CorrelationData, i :: Int64, j :: Int64, k :: Int64, l :: Int64)
  i = iseq_from_resnum(Correlation,i)
  j = iseq_from_resnum(Correlation,j)
  k = iseq_from_resnum(Correlation,k)
  l = iseq_from_resnum(Correlation,l)
  if min(i,j) == min(k,l) && max(i,j) == max(k,l)
    return 1., 1.
  else
    index = findi4(Correlation,i,j,k,l)
    if index == nothing
      return nothing, nothing
    else
      return Correlation.dist[index], Correlation.bin[index]
    end
  end
end

#
# Get the correlation of a constraint [i,j] with all other constraints defined by [k,X]
#
function getcorr(Correlation :: CorrelationData, i :: Int64, j :: Int64, k :: Int64)

  i = iseq_from_resnum(Correlation,i)
  j = iseq_from_resnum(Correlation,j)
  k = iseq_from_resnum(Correlation,k)

  nval = 0
  for l in 1:Correlation.nCA 
    val = getcorr(Correlation,i,j,k,l)
    if val != nothing
      nval = nval + 1
    end
  end
  index = Vector{Int64}(undef,nval)
  dist_corr = Vector{Float64}(undef,nval)
  bin_corr = similar(dist_corr)
  ival = 0
  for l in 1:Correlation.nCA
    dval, binval  = getcorr(Correlation,i,j,k,l)
    if dval != nothing
      ival = ival + 1
      index[ival] = l
      dist_corr[ival] = dval
      bin_corr[ival] = binval
    end
  end
  return index, dist_corr, bin_corr

end

#
# Get all the correlations of a set of constraints
#

function getcorr(Correlation :: CorrelationData, set :: Matrix{Int64}; self = false)

  nset = size(set)[1]
  npairs = Int64(nset*(nset-1)/2)
  if ! self
    iself = 1
  else
    iself = 0
    npairs = npairs + nset
  end
  index = Matrix{Int64}(undef,npairs,4)
  dist_corr = Vector{Float64}(undef,npairs)
  bin_corr = Vector{Float64}(undef,npairs)
  ipair = 0
  for i in 1:nset-iself
    i1 = iseq_from_resnum(Correlation,set[i,1])
    i2 = iseq_from_resnum(Correlation,set[i,2])
    for j in i+iself:nset
      ipair = ipair + 1
      j1 = iseq_from_resnum(Correlation,set[j,1])
      j2 = iseq_from_resnum(Correlation,set[j,2])
      dval, binval = getcorr(Correlation,i1,i2,j1,j2)
      if dval == nothing
        error(" Could not find correlation for: [",set[i,1],",",set[i,2],"]-[",set[j,1],",",set[j,2],"]")
      end
      index[ipair,1] = i1
      index[ipair,2] = i2
      index[ipair,3] = j1
      index[ipair,4] = j2
      dist_corr[ipair] = dval 
      bin_corr[ipair] = binval 
    end
  end
  return index, dist_corr, bin_corr

end

#
# Get the correlations between all pairs of two sets of constraints
#

function getcorr(Correlation :: CorrelationData, set1 :: Matrix{Int64}, set2 :: Matrix{Int64})

  nset1 = size(set1)[1]
  nset2 = size(set2)[2]
  npairs = nset1*nset2
  index = Matrix{Int64}(undef,npairs,4)
  dist_corr = Vector{Float64}(undef,npairs)
  bin_corr = Vector{Float64}(undef,npairs)
  ipair = 0
  for i in 1:nset1
    i1 = iseq_from_resnum(Correlation,set1[i,1])
    i2 = iseq_from_resnum(Correlation,set1[i,2])
    for j in 1:nset2
      ipair = ipair + 1
      j1 = iseq_from_resnum(Correlation,set2[j,1])
      j2 = iseq_from_resnum(Correlation,set2[j,2])
      dval, binval = getcorr(Correlation,i1,i2,j1,j2)
      if dval == nothing
        error(" Could not find correlation for: [",set1[i,1],",",set1[i,2],"]-[",set2[j,1],",",set2[j,2],"]")
      end
      index[ipair,1] = i1
      index[ipair,2] = i2
      index[ipair,3] = j1
      index[ipair,4] = j2
      dist_corr[ipair] = dval 
      bin_corr[ipair] = binval 
    end
  end
  return index, dist_corr, bin_corr

end










