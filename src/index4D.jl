
function index4D(dim,i,j,k,l)

  i1 = min(i,j)
  i2 = max(i,j)
  j1 = min(k,l)
  j2 = max(k,l)
  if ( i1 < j1 ) || ( i1 == j1 && i2 < j2 ) 
    i = i1
    j = i2
    k = j1
    l = j2
  else
    i = j1
    j = j2
    k = i1
    l = i2
  end

  dim2 = dim^2
  dim3 = dim2*dim
  index = (l-1)*dim3 + (k-1)*dim2 + (j-1)*dim + i
  return index

end
