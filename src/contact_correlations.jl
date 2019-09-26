
#Requires:
#using Statistics
#using JPleXL # from https://github.com/mcubeg/JPleXL
#using ProgressMeter

d(x,y) = sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2 )

function contact_correlations(pdblist :: Vector{String}; 
                              lastpdb=nothing, minsep=1, reference=nothing, tol=2.0)

  if lastpdb == nothing
    npdbs = length(pdblist)
  else
    npdbs = lastpdb
  end
  println(" Number of PDBs: ", npdbs)
  println(" Reading PDBs and computing contacts ... ")

  # If no external reference structure was provided, use the first frame
  if reference == nothing
    println(" Binary correlations will be computed using reference: ", pdblist[1])
    reference = pdblist[1]
  end
  println(" Binary tolerance = ", tol)
  pdb_ref = readPDB(reference)
  cas_ref = JPleXL.xCA(pdb_ref)
  ncas = size(cas_ref)[1]
  println(" Number of CAs: ", ncas)

  # Vector that contains the residue number for each CA in sequence: 
  resnum_from_seq = zeros(Int64,ncas)
  iseq = 0
  for iatom in 1:length(pdb_ref)
    if pdb_ref[iatom].name == "CA"
      iseq = iseq + 1 
      resnum_from_seq[iseq] = pdb_ref[iatom].resnum
    end
  end

  # If no contact list was provided, compute all contacts
  if minsep < 1 
    println(" ERROR: minsep must be greater or equal than 1 ")
  end
  ncontacts = Int64(ncas*(ncas-1)/2)
  for i in 1:minsep-1
    ncontacts = ncontacts - (ncas-i)
  end
 
  # Reference distances
  dref = Vector{Float64}(undef,ncontacts)

  # Setting array of contact indexes
  Contacts = Matrix{Int64}(undef,ncontacts,2)
  icontact = 0
  for i in 1:ncas-minsep
    for j in i+minsep:ncas
      icontact = icontact + 1
      Contacts[icontact,1] = i
      Contacts[icontact,2] = j
      dref[icontact] = d(cas_ref[i,1:3],cas_ref[j,1:3])
    end
  end
  println(" Number of contacts (i-j>",minsep,"): ", ncontacts)
  
  # Computing contact distances for all PDBs
  println(" Computing contact distances for all PDBs ... ")
  ContactDistances = Matrix{Float64}(undef,npdbs,ncontacts)
  ContactBin = Matrix{Bool}(undef,npdbs,ncontacts)
  @showprogress for ipdb in 1:npdbs
    pdb = readPDB(pdblist[ipdb])
    cas = JPleXL.xCA(pdb)
    for icontact in 1:ncontacts
      i = Contacts[icontact,1]
      j = Contacts[icontact,2]
      ContactDistances[ipdb,icontact] = d(cas[i,1:3],cas[j,1:3])
      if abs(ContactDistances[ipdb,icontact] - dref[icontact]) < tol
        ContactBin[ipdb,icontact] = true
      else
        ContactBin[ipdb,icontact] = false
      end
    end
  end
  for i in 1:ncas
    if resnum_from_seq[i] == 0
      error(" A residue in the sequence has no resnum associated with it: ", i)
    end
  end

  # Computing correlations
  println(" Computing correlations between contacts ... ")
  ncontact_pairs = Int64(ncontacts*(ncontacts-1)/2)
  index = Vector{Int64}(undef,ncontact_pairs)
  dist = Vector{Float64}(undef,ncontact_pairs)
  bin = Vector{Float64}(undef,ncontact_pairs)
  ipair = 0
  @showprogress for i in 1:ncontacts-1
    i1 = Contacts[i,1] 
    i2 = Contacts[i,2]
    for j in i+1:ncontacts
      j1 = Contacts[j,1]
      j2 = Contacts[j,2]
      ipair = ipair + 1
      index[ipair] = index4D(ncas,i1,i2,j1,j2) 
      dist[ipair] = Statistics.cor(ContactDistances[:,i],ContactDistances[:,j])
      # The binary correlation can only be computed with the binary count
      # varies within models (if there is variance)
      n_i_true = count(ContactBin[:,i])
      n_j_true = count(ContactBin[:,j])
      var_i = Statistics.std(ContactBin[:,i])
      var_j = Statistics.std(ContactBin[:,j])
      if var_i == 0. || var_j == 0.
        bin[ipair] = 0.
      else
        cov = Statistics.cov(ContactBin[:,i],ContactBin[:,j])
        bin[ipair] = cov / ( var_i * var_j )
      end
    end
  end

  println(" Sorting correlations by 4D matrix indexes... ")
  order = Vector{Int64}(undef,ncontact_pairs)
  for i in 1:ncontact_pairs
    order[i] = i
  end
  order = sort!(order,by=i->index[i])
  index = index[order]
  dist = dist[order]
  bin = bin[order]

  Correlation = CorrelationData(ncas,
                                ncontact_pairs,
                                ncontacts,
                                index,
                                dist,
                                bin,
                                resnum_from_seq,
                                ContactDistances,
                                ContactBin,
                                Contacts,
                                dref)

  return Correlation

end                            


