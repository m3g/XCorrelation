#
# Requires:
#
# Statistics
# JPleXL # from https://github.com/mcubeg/JPleXL
# ProgressMeter

using XCorrelation ; const XC = XCorrelation

# PDB file:
pdb="/home/guilherme/DCA_Modeling/SALB3.NOCST.001/pdb/S_00000001.pdb"

# Compute information
#    result: a square NxN matrix with the information provided by each contact
# 

information = XC.censoni_information(pdb)



