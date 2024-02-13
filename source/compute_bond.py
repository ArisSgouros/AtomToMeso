import os
import sys
import ast
import numpy as np
import math as m

DUMP_COL_ID    = 0
DUMP_COL_MOLID = 1
DUMP_COL_TYPE  = 2
DUMP_COL_X     = 3
DUMP_COL_Y     = 4
DUMP_COL_Z     = 5

if len(sys.argv) < 2:
   print ( )
   print ( "The required formats are the following:" )
   print ( )
   print ( "python lmp_acf.py \"LMP_DATAFILE\" \"[BondTypes1, BondType2, ..]\"" )
   print ( )
   print ( "example: python lmp_bond_dist.py LMP_DATAFILE [1,3]" )
   print ( "(this calculates the connectivity of bonds with type 1 and 3" )
   print ( "and outputs the atom IDs to be dumped from MD simulations in LAMMPS)" )
   print ( )
   print ( "or" )
   print ( )
   print ( "python lmp_acf.py \"LMP_DATAFILE\" \"[BondTypes1, BondType2, ..]\" \"NumberOfFrames\" \"LMP_TRJ_FILE\" " )
   print ( "example: python lmp_bond_dist.py LMP_DATAFILE [3] 21 LMP_TRJ_FILE" )
   print ( "(this estimated the bond length distribution)" )
   print ( )
   print ( "exiting.." )
   sys.exit()


#Set bond type(s)
btypes = ast.literal_eval(sys.argv[2])
btypes = [int(ii) for ii in btypes]

print ( "Selected bond types:",btypes )
FDESCRIPTION = ("_").join([ str(i) for i in btypes])

fileLines = []
atomIDs = []
bonds = []

#
# Read the DATA file and get the connectivity
#
print ( "Reading the connectivity from:",sys.argv[1],".." )

iLine = 0
with open(sys.argv[1],"r") as openFileObject:
   for curLine in openFileObject:
      if "bonds" in curLine:
         nBonds = int(curLine.split()[0])
      if "Bonds" in curLine:
         bondLineStart = iLine + 2
         bondLineEnd = bondLineStart + nBonds

      fileLines.append(curLine)
      iLine += 1

for curLine in range(bondLineStart,bondLineEnd):
   curLineSplit = fileLines[curLine].split()
   if int(curLineSplit[1]) in btypes:
      bonds.append([int(curLineSplit[2]),int(curLineSplit[3])])
      atomIDs.append(int(curLineSplit[2]))
      atomIDs.append(int(curLineSplit[3]))

#print ( bonds )

nBond = len(bonds)
print ( "A total of",nBond,"bonds were selected" )
#
# Print the necessary atom IDs for LAMMPS
#
#Sort the list and remove duplicates if any
atomIDs = sorted(set(atomIDs))
#Write the atom IDs in a file
f = open('o.'+FDESCRIPTION+'.atomIDs.dat', 'w')
for ID in atomIDs:
    f.write(str(ID)+" ")
f.close()
print ( "Printing the atom IDs in o.atomIDs.dat.." )
#
# Print the connectivity
#
f = open('o.'+FDESCRIPTION+'.bondIDs.dat', 'w')
for bond in bonds:
    f.write(str(bond[0]) + " " + str(bond[1]) + "\n")
f.close()
print ( "Printing the connectivity in o.bondIDs.dat.." )
#
# In case the configuration file is supplied proceed
# on the calculation of bonds
#
if len(sys.argv) == 3:
    sys.exit(0)

#Set the number of Frames
nFrame = int(sys.argv[3])
print ( "Number of Frames:",nFrame )

# Initialize the array of vectors with dimensions:
# [Nframe x NBonds x 3]
segvec = [[[0.0 for d in range(3)] for j in range(nBond)] for t in range(nFrame)]
#
# Load the atom trajectories
#
f = open(sys.argv[4],"r")
print ( "Reading the trajectory file",sys.argv[4],".." )
#for NFrame in range(20):
for tt in range(nFrame):

   f.readline()                # ITEM: TIMESTEP
   Timestep = int(f.readline())

   if (nFrame > 10 and tt % int(nFrame / 10.0) == 0):
      print ( "time step = ", Timestep )

   f.readline()                # ITEM: NUMBER OF ATOMS
   nAtoms = int(f.readline())
   f.readline()                # ITEM: BOX BOUNDS
   Dim = f.readline().split()  # xlo xhi
   LL0 = float(Dim[1])-float(Dim[0])
   Dim = f.readline().split()  # ylo yhi
   LL1 = float(Dim[1])-float(Dim[0])
   Dim = f.readline().split()  # zlo zhi
   LL2 = float(Dim[1])-float(Dim[0])
   f.readline()                # ITEM: FORMAT
   #
   # Loop over all atoms and create a
   # dictionary for atom IDs and coordinates
   #
   pos = {}
   for ii in range(nAtoms):
      line = f.readline().split()
      id = int(line[DUMP_COL_ID])
      xx = float(line[DUMP_COL_X])
      yy = float(line[DUMP_COL_Y])
      zz = float(line[DUMP_COL_Z])
      pos.update({id:[xx,yy,zz]})
   #
   # Now lets calculate the bond vectors!
   #
   bId = 0
   for bond in bonds:

       bondLen0 = pos[bond[1]][0]-pos[bond[0]][0]
       bondLen1 = pos[bond[1]][1]-pos[bond[0]][1]
       bondLen2 = pos[bond[1]][2]-pos[bond[0]][2]

       # Lets perform a minimum image just to be sure!
       bondLen0 -= LL0 * round(bondLen0 / LL0)
       bondLen1 -= LL1 * round(bondLen1 / LL1)
       bondLen2 -= LL2 * round(bondLen2 / LL2)

       segvec[tt][bId]=[bondLen0,bondLen1,bondLen2]

       bId += 1
f.close()
#
# Get the bond length distributions
#

seg_len = [[None for j in range(nBond)] for j in range(nFrame)]
#seg_len = []
for bId in range(nBond):
   for tt in range(nFrame):
      bond_len = m.sqrt(   segvec[tt][bId][0] * segvec[tt][bId][0] \
                         + segvec[tt][bId][1] * segvec[tt][bId][1] \
                         + segvec[tt][bId][2] * segvec[tt][bId][2] )

      seg_len[tt][bId] = bond_len

all_segs = [ item for sublist in seg_len for item in sublist ]

lbin = 0.006525
bins=[lbin*ii for ii in range(int(max(all_segs)/lbin)+1)]
nbins = len(bins)-1
st = np.histogram(all_segs, bins=bins)
norm_factor = lbin * np.sum(st[0])

f = open("o."+FDESCRIPTION+".bond_dist.dat","w")
f.write( " AVE: %16.9f \n" % np.average(all_segs))
f.write( " STD: %16.9f \n" % np.std(all_segs))
for ibin in range(0,nbins):
   f.write( "%16.9f  %16.9f \n" % (st[1][ibin]+0.5*lbin, (float(st[0][ibin]))/norm_factor))
   #f.write( "%16.9f  %16.9f \n" % (st[1][ibin], (float(st[0][ibin]))/norm_factor))
f.close()

#
# Get the partial bond length distributions
#

st = [None] * nBond

for aid in range(nBond):
   aux = []
   for tt in range(nFrame):
      aux.append(seg_len[tt][aid])
   st[aid] = np.histogram(aux, bins=bins)

f = open("o."+FDESCRIPTION+".bond_dist_partial.dat","w")
for ibin in range(nbins):
   f.write( "%16.9f " % (st[0][1][ibin]+0.5*lbin))
   for aid in range(nBond):
      f.write( "%d " % (st[aid][0][ibin]))
   f.write( "\n" )
f.close()
