#
# Requires:
#
# using Statistics
# using JPleXL # from https://github.com/mcubeg/JPleXL
# using ProgressMeter

push!(LOAD_PATH,"/home/leandro/programs/JPleXL") # from https://github.com/mcubeg/JPleXL
push!(LOAD_PATH,"/home/leandro/Drive/Alunos/Guilherme/wdegree/XCorrelation") 

using XCorrelation ; const XC = XCorrelation

# Read list of PDB files
pdbdir = "/home/guilherme/DCA_Modeling/SALB3.NOCST.001/pdb"
println(" PDB Directory: ", pdbdir)
pdblist = readdir(pdbdir)
pdblist = filter( x -> occursin("S_",x),pdblist)
@. pdblist = pdbdir*"/"*pdblist

# Reference PDB file:
reference="/home/guilherme/DCA_Modeling/SALB3.NOCST.001/pdb/S_00000001.pdb"

# Compute correlations
#    use: minsep=1 to compute for all pairs of residues (i-j>1) (default if argument is removed)
#         lastpdb=nothing to compute for all PDBs (default if argument is removed)

C = XC.contact_correlations(pdblist,lastpdb=5,tol=1.0,reference=reference,minsep=100)

#
# Results: 
#

# C.Contacts contains the indexes of the residues involved in each contact pair

# C.ContactDref contins the reference distances for each pair

# C.ContactBin is a Ncontact x Npdb matrix containing 0 or 1 for the satisfaction of the
#              contact at each PDB file (with tolerance +/- tol)

# C.ContactDistances is a Ncontact x Npdb matrix containing the correlation of the distances
#                    of the contacts


