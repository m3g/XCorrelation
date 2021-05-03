using SpecialFunctions


euclidean_distance(x, y) = sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2 )
gamma_inc_prec(x, y) = gamma(x) * (1.0-gamma_inc(x, y, 0)[1]) # high precision, first output, matches Wolfram definition


function F_censoni(r, l)
  # primitive
  Al2v = 16.58 * (l^1.15)
  B = 2.4343
  C = 0.8555
  D = -0.3125
  E = 0.9375
  return 1.0 + (B * D * gamma_inc_prec( E, (C/Al2v)^1.6 * r^3.2)) / C^1.5
end


function f_censoni(rmin, rmax, l)
  return -log2( F_censoni( rmax, l) - F_censoni(rmin, l))
end


function f_censoni(rmin, rmax, l, stp)
  # numerical integration, reasonable step is: stp = 0.0001
  Al2v = 16.58 * (l^1.15)
  f = 0.
  s = rmin
  while s < rmax
    r2 = s^2
    f = f +( 2.4343 / Al2v^(1.5) )*exp( -(0.8555 * r2/Al2v)^1.6) * r2 * stp
    s = s + stp
  end
  return -log2(f)
end


function censoni_information(pdbfile; select = "name CA")
  pdb = PDBTools.readPDB(pdbfile)
  cas = PDBTools.coor(pdb, select; xyz_in_cols = true)

  # example use case, but should be changed
  # to receive confidence interval as input
  precision = 0.5

  nCA = size(cas)[1]
  information = zeros(nCA,nCA)

  for i in 1:nCA-1
    for j in i+1:nCA
      r = euclidean_distance(cas[i,:],cas[j,:])
      l = j-i
      information[i,j] = f_censoni( r-precision, r+precision, l)
      information[j,i] = information[i,j]
    end
  end

  return information
end
