#
# Computes the mean correlation of a set of constraints
#
function meancorr(Correlation :: CorrelationData, set :: Matrix{Int64})
  nconst = size(set)[1]
  av_corr = 0.
  for i in 1:nconst-1
    i1 = set[i,1]
    i2 = set[i,2]
    for j in i+1:nconst
      j1 = set[j,1]
      j2 = set[j,2]
      av_corr = av_corr + getcorr(Correlation,i1,i2,j1,j2)
    end
  end
  av_corr = av_corr / ( nconst*(nconst-1)/2 )
  return av_corr
end
