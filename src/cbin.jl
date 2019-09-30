# Function to compute the binary correlations from the already computed distance data 

function cbin(C :: CorrelationData; tol=2.0)
  
  npdbs = size(ContactBin)[1]
  Cbin = Matrix{Int64}(undef,npdbs,C.ncontacts)
  for ipdb in 1:npdbs
    for ic in 1:C.ncontacts
      if abs(C.ContactDistances[ipdb,ic] - C.ContactDref[ic]) <= tol
        Cbin[ipdb,ic] = 1
      else
        Cbin[ipdb,ic] = 0
      end
    end
  end
  return Cbin

end







