import ast
import numpy as np
import math as m
from export import GetDumpFormat

def ComputeAngle(path_data, types, nFrame, path_dump, export_hist=False, export_hist_partial=False, lbin=0.1, verbose=False, debug=False):

   # Parse the format of the dump file
   dump_col = GetDumpFormat(path_dump)
   DUMP_COL_ID = dump_col['id']
   DUMP_COL_X = dump_col['x']
   DUMP_COL_Y = dump_col['y']
   DUMP_COL_Z = dump_col['z']

   if verbose: print( "Selected types:",types )
   fname = ("_").join([ str(i) for i in types])

   lines = []
   atoms = []
   angles = []

   #
   # Read the DATA file and get the connectivity
   #
   if verbose: print( "Reading the connectivity from:",path_data,".." )

   iline = 0
   with open(path_data,"r") as foo:
      for cur_line in foo:
         if "angles" in cur_line:
            n_angle_tot = int(cur_line.split()[0])
         if "Angles" in cur_line:
            angle_line_start = iline + 2
            angle_line_end = angle_line_start + n_angle_tot

         lines.append(cur_line)
         iline += 1

   for cur_line in range(angle_line_start,angle_line_end):
      lsplit = lines[cur_line].split()
      if int(lsplit[1]) in types:
         angles.append([int(lsplit[2]),int(lsplit[3]),int(lsplit[4])])
         atoms.append(int(lsplit[2]))
         atoms.append(int(lsplit[3]))
         atoms.append(int(lsplit[4]))

   n_count = len(angles)
   if verbose: print( "A total of",n_count,"angles were selected" )
   #
   # Print the necessary atom IDs for LAMMPS
   #
   #Sort the list and remove duplicates if any
   atoms = sorted(set(atoms))
   
   if debug:
      print( "Printing the atom IDs in o.atomIDs.dat.." )
      f = open('o.'+fname+'.atomIDs.dat', 'w')
      for ID in atoms:
         f.write(str(ID)+" ")
      f.close()
      print( "Printing the connectivity in o.angleIDs.dat.." )
      f = open('o.'+fname+'.angleIDs.dat', 'w')
      for angle in angles:
         f.write(str(angle[0]) + " " + str(angle[1]) + " " + str(angle[2]) + "\n")
      f.close()

   if verbose: print( "Number of Frames:",nFrame )

   # Initialize the array of vectors with dimensions:
   # [Nframe x NAngles x 3]
   val_indiv_tt = [[0.0 for j in range(n_count)] for t in range(nFrame)]
   val_all = []
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
      for angle in angles:

          ra0 = pos[angle[0]][0]-pos[angle[1]][0]
          ra1 = pos[angle[0]][1]-pos[angle[1]][1]
          ra2 = pos[angle[0]][2]-pos[angle[1]][2]

          # Lets perform a minimum image just to be sure!
          ra0 -= LL0 * round(ra0 / LL0)
          ra1 -= LL1 * round(ra1 / LL1)
          ra2 -= LL2 * round(ra2 / LL2)

          ra_mag = m.sqrt(ra0*ra0+ra1*ra1+ra2*ra2)

          rb0 = pos[angle[2]][0]-pos[angle[1]][0]
          rb1 = pos[angle[2]][1]-pos[angle[1]][1]
          rb2 = pos[angle[2]][2]-pos[angle[1]][2]

          # Lets perform a minimum image just to be sure!
          rb0 -= LL0 * round(rb0 / LL0)
          rb1 -= LL1 * round(rb1 / LL1)
          rb2 -= LL2 * round(rb2 / LL2)

          rb_mag = m.sqrt(rb0*rb0+rb1*rb1+rb2*rb2)

          val_indiv_tt[tt][aid] = m.acos( (ra0*rb0 + ra1*rb1+ ra2*rb2)  / (ra_mag * rb_mag) )
          aid += 1
   f.close()

   val_all = [ item for sublist in val_indiv_tt for item in sublist ]

   mean = np.average(val_all)
   std = np.std(val_all)

   if export_hist:
      bins=[lbin*ii for ii in range(0, int(max(val_all)/lbin)+1)]
      nbins = len(bins)-1
      st = np.histogram(val_all, bins=bins)
      norm_factor = lbin * np.sum(st[0])

      f = open("o."+fname+".angle_dist.dat","w")
      f.write( " AVE: %16.9f \n" % np.average(val_all))
      f.write( " STD: %16.9f \n" % np.std(val_all))
      for ibin in range(nbins):
         f.write( "%16.9f  %16.9f \n" % (st[1][ibin]+0.5*lbin, float(st[0][ibin])/norm_factor))
      f.close()

   if export_hist_partial:
      st = [None] * n_count
      for aid in range(n_count):
         aux = []
         for tt in range(nFrame):
            aux.append(val_indiv_tt[tt][aid])
         st[aid] = np.histogram(aux, bins=bins)

      f = open("o."+fname+".angle_dist_partial.dat","w")
      for ibin in range(nbins):
         f.write( "%16.9f " % (st[0][1][ibin]+0.5*lbin))
         for aid in range(n_count):
            f.write( "%d " % (st[aid][0][ibin]))
         f.write( "\n" )
      f.close()

   return mean, std
