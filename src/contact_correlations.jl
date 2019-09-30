
function contact_correlations(pdblist :: Vector{String}; 
                              lastpdb=nothing, minsep=1, reference=nothing, tol=2.0)

  C = contact_data(pdblist, lastpdb=lastpdb, minsep=minsep, reference=reference, tol=tol)

  contatct_correlation!(C)
  
  return C

end

function contact_correlation!(C :: CorrelationData)

  # Computing correlations

  println(" Computing correlations between contacts ... ")

  ipair = 0
  @showprogress for i in 1:C.ncontacts-1
    i1 = C.Contacts[i,1] 
    i2 = C.CContacts[i,2]
    for j in i+1:C.ncontacts
      j1 = C.Contacts[j,1]
      j2 = C.Contacts[j,2]
      ipair = ipair + 1
      C.index[ipair] = index4D(C.nCA,i1,i2,j1,j2) 
      C.dist[ipair] = Statistics.cor(C.ContactDistances[:,i],C.ContactDistances[:,j])
      # The binary correlation can only be computed with the binary count
      # varies within models (if there is variance)
      n_i_true = count(C.ContactBin[:,i])
      n_j_true = count(C.ContactBin[:,j])
      var_i = Statistics.std(C.ContactBin[:,i])
      var_j = Statistics.std(C.ContactBin[:,j])
      if var_i == 0. || var_j == 0.
        C.bin[ipair] = 0.
      else
        cov = Statistics.cov(C.ContactBin[:,i],C.ContactBin[:,j])
        C.bin[ipair] = cov / ( var_i * var_j )
      end
    end
  end

  println(" Sorting correlations by 4D matrix indexes... ")
  order = Vector{Int64}(undef,C.npairs)
  for i in 1:C.npairs
    order[i] = i
  end
  order = sort!(order,by=i->C.index[i])
  C.index = C.index[order]
  C.dist = C.dist[order]
  C.bin = C.bin[order]

  return C

end                            



