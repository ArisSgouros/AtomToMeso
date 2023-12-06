from copy import deepcopy
from network import EMPTY, c_atom

def ExportLammpsData(filename, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box):
   f = open(filename, 'w')

   n_beads = 0
   for group in network:
      n_beads += len(group.beads)

   # get the number of bead types
   aux = set()
   for type in num_of_type.values():
      aux.add(type)
   n_bead_types = len(aux)

   f.write('# Coarse grained LAMMPS data file\n')
   f.write('\n')
   f.write('%d atoms\n' % n_beads)
   f.write('%d bonds\n' % len(bead_bonds))
   f.write('%d angless\n' % len(bead_angles))
   f.write('\n')
   f.write('%d atom types\n' % len(mass_of_bead))
   f.write('%d bond types\n' % len(num_of_bond_type))
   f.write('%d angle types\n' % len(num_of_angle_type))
   f.write('\n')
   f.write('%s %s xlo xhi\n' % (box['xlo'], box['xhi']))
   f.write('%s %s ylo yhi\n' % (box['ylo'], box['yhi']))
   f.write('%s %s zlo zhi\n' % (box['zlo'], box['zhi']))
   f.write('\n')
   f.write('Atoms\n')
   f.write('\n')

   for group in network:
      for bead in group.beads.values():
         f.write("%d %d %d %16.9f %16.9f %16.9f %16.9f # %s\n" % ( bead.bId, bead.molId, num_of_type[bead.type], bead.q, bead.x, bead.y, bead.z, bead.type ))

   if bead_bonds:
      f.write('\n')
      f.write('Bonds\n')
      f.write('\n')
      for bond in bead_bonds:
         f.write("%d %d %d %d # %s_%s\n" % ( bond.Id, bond.type, bond.i, bond.j, bond.itype, bond.jtype ))

   if bead_angles:
      f.write('\n')
      f.write('Angles\n')
      f.write('\n')
      for angle in bead_angles:
         f.write("%d %d %d %d %d # %s_%s_%s\n" % ( angle.Id, angle.type, angle.i, angle.j, angle.k, angle.itype, angle.jtype, angle.ktype ))

   f.write('\n')
   f.write('Masses\n')
   f.write('\n')
   for bead_type in mass_of_bead:
      f.write("%d %16.9f # %s\n" % (num_of_type[bead_type], mass_of_bead[bead_type], bead_type))

   f.write('\n')
   f.write('Pair Coeffs\n')
   f.write('\n')
   for bead_type in mass_of_bead:
      f.write("%s DPD_COEFF_OF_%s # %s\n" % (num_of_type[bead_type], bead_type, bead_type))

   if num_of_bond_type:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for bond in num_of_bond_type:
         f.write(str(num_of_bond_type[bond])+" Bond coeffs # " + bond+"\n")

   if num_of_angle_type:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for angle in num_of_angle_type:
         f.write(str(num_of_angle_type[angle])+" Angle coeffs # " + angle+"\n")

   f.close()



def ExportLammpsDump(path_dump_in, path_dump_out, N_FRAME, EV_FRAME, network, mass_of_type, dump_col, num_of_type):

   DUMP_COL_ID = dump_col['id']
   DUMP_COL_MOLID = dump_col['molid']
   DUMP_COL_TYPE = dump_col['type']
   DUMP_COL_X = dump_col['x']
   DUMP_COL_Y = dump_col['y']
   DUMP_COL_Z = dump_col['z']


   fin = open(path_dump_in, 'r')
   fout = open(path_dump_out, 'w')

   n_beads = 0
   for group in network:
      n_beads += len(group.beads)

   iframe = -1

   while True:

      iframe += 1

      if N_FRAME != "MAX":
         if (iframe > 10 and iframe % int(N_FRAME / 10.0) == 0):
            print(iframe)
         if iframe == N_FRAME:
            break

      if (iframe % EV_FRAME != 0):
         continue

      line = fin.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break
      fout.write(line)

      line = fin.readline()
      fout.write(line)

      line = fin.readline() # ITEM: NUMBER OF ATOMS
      fout.write(line)

      n_atoms = int(fin.readline().split()[0]) #

      fout.write(str(n_beads)+"\n")
      line = fin.readline() # ITEM: BOX BOUNDS pp pp pp
      fout.write(line)
      line = fin.readline() # X
      fout.write(line)
      line = fin.readline() # Y
      fout.write(line)
      line = fin.readline() # Z
      fout.write(line)
      line = fin.readline() # ITEM: ATOMS id mol type xu yu zu
      fout.write(line)

      atom_list = {}
      for ii in range(n_atoms):
         line_split = fin.readline().split()
         iat = c_atom()
         iat.Id    = int(  line_split[DUMP_COL_ID])
         iat.molId = int(  line_split[DUMP_COL_MOLID])
         iat.type  = int(  line_split[DUMP_COL_TYPE])
         iat.x     = float(line_split[DUMP_COL_X])
         iat.y     = float(line_split[DUMP_COL_Y])
         iat.z     = float(line_split[DUMP_COL_Z])
         atom_list[iat.Id] = deepcopy(iat)


      for group in network:
         for bead in group.beads.values():
            bead_x = 0
            bead_y = 0
            bead_z = 0
            bead_mass = 0
            for aId in bead.aIds:
               atom = atom_list[aId]
               imass = mass_of_type[atom.type]

               bead_x += atom.x * imass
               bead_y += atom.y * imass
               bead_z += atom.z * imass
               bead_mass += imass

            bead.x = bead_x / bead_mass
            bead.y = bead_y / bead_mass
            bead.z = bead_z / bead_mass
            bead.mass = bead_mass

            fout.write("%d %d %d %16.9f %16.9f %16.9f\n" %(bead.bId, bead.molId, num_of_type[bead.type], bead.x, bead.y, bead.z))


   fin.close()
   fout.close()
