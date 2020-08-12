
d(x,y) = sqrt( (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2 )

function contact_data(pdblist :: Vector{String}; 
                      lastpdb=nothing, minsep=1, reference=nothing, tol=2.0, correlations=false,
                      selection="name CA")

  if lastpdb == nothing
    npdbs = length(pdblist)
  else
    npdbs = lastpdb
  end
  println(" Number of PDBs: ", npdbs)
  println(" Reading PDBs and computing contacts ... ")

  # If no external reference structure was provided, use the first frame

  if typeof(reference) == Float64
    println(" Contacts will be computed if distance is lower than: ", reference)
    ref_dist = reference
    reference = pdblist[1]
  elseif reference == nothing
    println(" Binary correlations will be computed using reference: ", pdblist[1])
    reference = pdblist[1]
  end
    println(" Binary tolerance = ", tol)
    pdb_ref = PDBTools.readPDB(reference)
    cas_ref = PDBTools.coor(pdb_ref, "name CA")
    nCA = size(cas_ref)[1]
    println(" Number of CAs: ", nCA)

  # Vector that contains the residue number for each CA in sequence: 
  resnum_from_seq = zeros(Int64,nCA)
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
  ncontacts = Int64(nCA*(nCA-1)/2)
  for i in 1:minsep-1
    ncontacts = ncontacts - (nCA-i)
  end
 
  # Reference distances
  Dref = Vector{Float64}(undef,ncontacts)

  # Setting array of contact indexes
  Contacts = Matrix{Int64}(undef,ncontacts,2)
  icontact = 0
  for i in 1:nCA-minsep
    for j in i+minsep:nCA
      icontact = icontact + 1
      Contacts[icontact,1] = i
      Contacts[icontact,2] = j
      Dref[icontact] = d(cas_ref[i,1:3],cas_ref[j,1:3])
    end
  end
  println(" Number of contacts (i-j>",minsep,"): ", ncontacts)
  
  # Computing contact distances for all PDBs
  println(" Computing contact distances for all PDBs ... ")
  ContactDistances = Matrix{Float64}(undef,npdbs,ncontacts)
  ContactBin = Matrix{Bool}(undef,npdbs,ncontacts)
  @showprogress for ipdb in 1:npdbs
    pdb = PDBTools.readPDB(pdblist[ipdb])
    cas = PDBTools.coor(pdb, selection)
    for icontact in 1:ncontacts
      i = Contacts[icontact,1]
      j = Contacts[icontact,2]
      ContactDistances[ipdb,icontact] = d(cas[i,1:3],cas[j,1:3])

    # Computing binary matrices based on a given reference PDB...
      if @isdefined(ref_dist)
        if (ContactDistances[ipdb,icontact] < ref_dist + tol)
          ContactBin[ipdb,icontact] = true
        else
          ContactBin[ipdb,icontact] = false
        end
      else 
        if abs(ContactDistances[ipdb,icontact] - Dref[icontact]) < tol
          ContactBin[ipdb,icontact] = true
        else
          ContactBin[ipdb,icontact] = false
        end
      end
    end
  end
  for i in 1:nCA
    if resnum_from_seq[i] == 0
      error(" A residue in the sequence has no resnum associated with it: ", i)
    end
  end

  npairs = Int64(ncontacts*(ncontacts-1)/2)
  if correlations
    index = Vector{Int64}(undef,npairs)
    dist = Vector{Float64}(undef,npairs)
    bin = Vector{Float64}(undef,npairs)
  else
    index = Vector{Int64}(undef,1)
    dist = Vector{Float64}(undef,1)
    bin = Vector{Float64}(undef,1)
  end

  C = CorrelationData(npdbs,
                      nCA,
                      npairs,
                      ncontacts,
                      index,
                      dist,
                      bin,
                      resnum_from_seq,
                      ContactDistances,
                      ContactBin,
                      Contacts,
                      Dref)

  return C

end                            


