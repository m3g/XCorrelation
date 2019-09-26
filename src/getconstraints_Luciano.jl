#
# Get the constraints from Luciano's file, given the indexes
# of the CA atoms, and returns the set of constraints as the
# index of residues
#
function getconstraints_Luciano(setfile_name,pdb)
  setfile = open(setfile_name,"r")
  nrest = 0
  for line in eachline(setfile)
    data = split(line)
    if data[1] != "#" && data[1] != "ATOM1" 
      i = parse(Int64,data[1])
      j = parse(Int64,data[2])
      nrest = nrest + 1
    end
  end
  close(setfile)
  set = Matrix{Int64}(undef,nrest,2)
  setfile = open(setfile_name,"r")
  #rewind(setfile)
  irest = 0
  for line in eachline(setfile)
    data = split(line)
    if data[1] != "#" && data[1] != "ATOM1" 
      i = parse(Int64,data[1])
      j = parse(Int64,data[2])
      irest = irest + 1
      for atom in pdb
        if atom.index == i  
          set[irest,1] = atom.resnum
          continue
        end
        if atom.index == j  
          set[irest,2] = atom.resnum
          continue
        end
      end
    end
  end
  close(setfile)
  return set
end














