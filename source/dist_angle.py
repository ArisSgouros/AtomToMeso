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
   print( )
   print( "The required formats are the following:" )
   print( )
   print( "python lmp_angle_disp.py \"LMP_DATAFILE\" \"[AngleTypes1, AngleType2, ..]\"" )
   print( )
   print( "example: python lmp_acf.py LMP_DATAFILE [1,3]" )
   print( "(this calculates the connectivity of angles with type 1 and 3" )
   print( "and outputs the atom IDs to be dumped from MD simulations in LAMMPS)" )
   print( )
   print( "or" )
   print( )
   print( "python lmp_acf.py \"LMP_DATAFILE\" \"[AngleTypes1, AngleType2, ..]\" \"NumberOfFrames\" \"LMP_TRJ_FILE\" " )
   print( "example: python lmp_angle_disp.py LMP_DATAFILE [3] 21 LMP_TRJ_FILE" )
   print( "(this estimated the angle length distribution)" )
   print( )
   print( "exiting.." )
   sys.exit()

#Set angle type(s)
atypes = ast.literal_eval(sys.argv[2])
atypes = [int(ii) for ii in atypes]

print( "Selected angle types:",atypes )
FDESCRIPTION = ("_").join([ str(i) for i in atypes])

fileLines = []
atomIDs = []
angle_ids = []

#
# Read the DATA file and get the connectivity
#
print( "Reading the connectivity from:",sys.argv[1],".." )

iLine = 0
with open(sys.argv[1],"r") as openFileObject:
   for curLine in openFileObject:
      if "angles" in curLine:
         nAngles = int(curLine.split()[0])
      if "Angles" in curLine:
         angleLineStart = iLine + 2
         angleLineEnd = angleLineStart + nAngles

      fileLines.append(curLine)
      iLine += 1

for curLine in range(angleLineStart,angleLineEnd):
   curLineSplit = fileLines[curLine].split()
   if int(curLineSplit[1]) in atypes:
      angle_ids.append([int(curLineSplit[2]),int(curLineSplit[3]),int(curLineSplit[4])])
      atomIDs.append(int(curLineSplit[2]))
      atomIDs.append(int(curLineSplit[3]))
      atomIDs.append(int(curLineSplit[4]))

nAngle = len(angle_ids)
print( "A total of",nAngle,"angles were selected" )
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
print( "Printing the atom IDs in o.atomIDs.dat.." )
#
# Print the connectivity
#
f = open('o.'+FDESCRIPTION+'.angleIDs.dat', 'w')
for angle in angle_ids:
    f.write(str(angle[0]) + " " + str(angle[1]) + " " + str(angle[2]) + "\n")
f.close()
print( "Printing the connectivity in o.angleIDs.dat.." )
#
# In case the configuration file is supplied proceed
# on the calculation of angles
#
if len(sys.argv) == 3:
    sys.exit(0)

#Set the number of Frames
nFrame = int(sys.argv[3])
print( "Number of Frames:",nFrame )

# Initialize the array of vectors with dimensions:
# [Nframe x NAngles x 3]
angles_frame = [[0.0 for j in range(nAngle)] for j in range(nFrame)]
all_angles = []
#
# Load the atom trajectories
#
f = open(sys.argv[4],"r")
print( "Reading the trajectory file",sys.argv[4],".." )
#for NFrame in range(20):
for tt in range(nFrame):

   f.readline()                # ITEM: TIMESTEP
   Timestep = int(f.readline())

   if (nFrame > 10 and tt % int(nFrame / 10.0) == 0):
      print( "time step = ", Timestep )

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
   # Now lets calculate the angle vectors!
   #
   bId = 0
   aid = 0
   for angle in angle_ids:

       bond_a0 = pos[angle[0]][0]-pos[angle[1]][0]
       bond_a1 = pos[angle[0]][1]-pos[angle[1]][1]
       bond_a2 = pos[angle[0]][2]-pos[angle[1]][2]

       # Lets perform a minimum image just to be sure!
       bond_a0 -= LL0 * round(bond_a0 / LL0)
       bond_a1 -= LL1 * round(bond_a1 / LL1)
       bond_a2 -= LL2 * round(bond_a2 / LL2)

       bond_a_mag = m.sqrt(bond_a0*bond_a0+bond_a1*bond_a1+bond_a2*bond_a2)

       bond_b0 = pos[angle[2]][0]-pos[angle[1]][0]
       bond_b1 = pos[angle[2]][1]-pos[angle[1]][1]
       bond_b2 = pos[angle[2]][2]-pos[angle[1]][2]

       # Lets perform a minimum image just to be sure!
       bond_b0 -= LL0 * round(bond_b0 / LL0)
       bond_b1 -= LL1 * round(bond_b1 / LL1)
       bond_b2 -= LL2 * round(bond_b2 / LL2)

       bond_b_mag = m.sqrt(bond_b0*bond_b0+bond_b1*bond_b1+bond_b2*bond_b2)

       angles_frame[tt][aid] = m.acos( (bond_a0*bond_b0 + bond_a1*bond_b1+ bond_a2*bond_b2)  / (bond_a_mag * bond_b_mag) )

       #print((bond_a0, bond_a1, bond_a2, bond_a_mag, bond_b0, bond_b1, bond_b2, bond_b_mag, "#", angle, ) )

   #    all_angles.append( angles_frame[tt][aid] )

       aid += 1
f.close()

all_angles = [ item for sublist in angles_frame for item in sublist ]

#
# Get the total angle length distributions
#

lbin = 0.01

#bins=[lbin*ii for ii in range(int(min(all_angles)/lbin)-2, int(max(all_angles)/lbin)+2)]
bins=[lbin*ii for ii in range(0, int(max(all_angles)/lbin)+2)]
nbins = len(bins)-1
st = np.histogram(all_angles, bins=bins)

norm_factor = lbin * np.sum(st[0])

f = open("o."+FDESCRIPTION+".angle_dist.dat","w")
f.write( " AVE: %16.9f \n" % np.average(all_angles))
f.write( " STD: %16.9f \n" % np.std(all_angles))
for ibin in range(nbins):
   f.write( "%16.9f  %16.9f \n" % (st[1][ibin]+0.5*lbin, float(st[0][ibin]) / norm_factor))
f.close()

#
# Get the partial angle length distributions
#

st = [None] * nAngle

for aid in range(nAngle):
   aux = []
   for tt in range(nFrame):
      aux.append(angles_frame[tt][aid])
   st[aid] = np.histogram(aux, bins=bins)

f = open("o."+FDESCRIPTION+".angle_dist_partial.dat","w")
for ibin in range(nbins):
   f.write( "%16.9f " % (st[0][1][ibin]+0.5*lbin))
   for aid in range(nAngle):
      f.write( "%d " % (st[aid][0][ibin]))
   f.write( "\n" )
f.close()


print( "Done!" )
sys.exit(0)
