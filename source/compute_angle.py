import os
import sys
import ast
import numpy as np
import math as m

def ComputeAngle(path_data, atypes, nFrame, path_dump, lbin, verbose, debug):
   DUMP_COL_ID    = 0
   DUMP_COL_MOLID = 1
   DUMP_COL_TYPE  = 2
   DUMP_COL_X     = 3
   DUMP_COL_Y     = 4
   DUMP_COL_Z     = 5

   if verbose: print( "Selected angle types:",atypes )
   FDESCRIPTION = ("_").join([ str(i) for i in atypes])

   fileLines = []
   atomIDs = []
   angle_ids = []

   #
   # Read the DATA file and get the connectivity
   #
   if verbose: print( "Reading the connectivity from:",path_data,".." )

   iLine = 0
   with open(path_data,"r") as openFileObject:
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
   if verbose: print( "A total of",nAngle,"angles were selected" )
   #
   # Print the necessary atom IDs for LAMMPS
   #
   #Sort the list and remove duplicates if any
   atomIDs = sorted(set(atomIDs))
   
   if debug:
      print( "Printing the atom IDs in o.atomIDs.dat.." )
      f = open('o.'+FDESCRIPTION+'.atomIDs.dat', 'w')
      for ID in atomIDs:
         f.write(str(ID)+" ")
      f.close()
      print( "Printing the connectivity in o.angleIDs.dat.." )
      f = open('o.'+FDESCRIPTION+'.angleIDs.dat', 'w')
      for angle in angle_ids:
         f.write(str(angle[0]) + " " + str(angle[1]) + " " + str(angle[2]) + "\n")
      f.close()

   if verbose: print( "Number of Frames:",nFrame )

   # Initialize the array of vectors with dimensions:
   # [Nframe x NAngles x 3]
   angles_frame = [[0.0 for j in range(nAngle)] for j in range(nFrame)]
   all_angles = []
   #
   # Load the atom trajectories
   #
   f = open(path_dump,"r")
   if verbose: print( "Reading the trajectory file",path_dump,".." )
   for tt in range(nFrame):

      f.readline()                # ITEM: TIMESTEP
      Timestep = int(f.readline())

      if verbose:
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
          aid += 1
   f.close()

   all_angles = [ item for sublist in angles_frame for item in sublist ]

   mean = np.average(all_angles)
   std = np.std(all_angles)
   #
   # Get the total angle length distributions
   #

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

   return mean, std
