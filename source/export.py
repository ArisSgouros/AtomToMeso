from copy import deepcopy
from network import Atom

def GetDumpFormat(filename):
   col = {}
   f = open(filename,"r")
   line = ""
   for ii in range(9):
      line = f.readline()
   f.close()
   header = line.split()
   header = header[2:] # pop first two elements ('ITEM:', 'ATOMS')
   for icol in range(len(header)):
      attrib = header[icol]
      if attrib == 'id':
         col['id'] = icol
      if attrib in ['type']:
         col['type'] = icol
      if attrib in ['mol']:
         col['molid'] = icol
      if attrib in ['x', 'xu', 'xs']:
         col['x'] = icol
      if attrib in ['y', 'yu', 'ys']:
         col['y'] = icol
      if attrib in ['z', 'zu', 'zs']:
         col['z'] = icol
   if not 'id' in col:
      print('col id not found in the dump file')
      sys.exit()
   if not 'x' in col:
      print('col x not found in the dump file')
      sys.exit()
   if not 'y' in col:
      print('col yy not found in the dump file')
      sys.exit()
   if not 'z' in col:
      print('col zz not found in the dump file')
      sys.exit()
   return col

def ParseLammpsData(path_lammps_data_in, network, atomtype):

   # Deal with the format of the Lammps data file
   if atomtype == 'full':
      data_col = {'id': 0, 'molid': 1, 'type': 2, 'q': 3, 'x':4, 'y':5, 'z':6, 'ix':7, 'iy':8, 'iz':9}
   else:
      print('unsupported format of Lammps data files for the atom type ' + atomtype)
      sys.exit()

   DATA_COL_ID    = data_col['id']
   DATA_COL_MOLID = data_col['molid']
   DATA_COL_TYPE  = data_col['type']
   DATA_COL_Q     = data_col['q']
   DATA_COL_X     = data_col['x']
   DATA_COL_Y     = data_col['y']
   DATA_COL_Z     = data_col['z']
   DATA_COL_IX    = data_col['ix']
   DATA_COL_IY    = data_col['iy']
   DATA_COL_IZ    = data_col['iz']

   gaId = 0
   for group in network:
     gaId += group.nat

   n_atoms = 0
   box = {}
   g = open(path_lammps_data_in, 'r')


   lines = []
   line_num = 0
   while True:
      line = g.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break

      line_split = line.split()

      if "atoms" in line:
         n_atoms = int(line_split[0])
         if n_atoms != gaId:
            print("ERROR: n_atoms (",n_atoms,") != gaId (",gaId,")")
            print("       there is a mismatch between the generated topology and the input data file")
            exit()

      if "atom types" in line:
         n_atom_types = int(line_split[0])

      if "xlo xhi" in line:
         box['xlo'] = line_split[0]
         box['xhi'] = line_split[1]
         Lx = float(box['xhi'])-float(box['xlo'])
      if "ylo yhi" in line:
         box['ylo'] = line_split[0]
         box['yhi'] = line_split[1]
         Ly = float(box['yhi'])-float(box['ylo'])
      if "zlo zhi" in line:
         box['zlo'] = line_split[0]
         box['zhi'] = line_split[1]
         Lz = float(box['zhi'])-float(box['zlo'])

      if "Atoms" in line:
         Atoms_start = line_num + 2

      if "Masses" in line:
         Masses_start = line_num + 2

      lines.append(line)
      line_num += 1

   g.close()

   mass_of_type = {}
   #
   # Read the Masses section
   print("Reading Masses section..")
   for ii in range(n_atom_types):
      line_split = lines[Masses_start + ii].split()
      type = int(line_split[0])
      mass = float(line_split[1])
      mass_of_type[type] = mass
      print("   ",type,mass)

   atom_list = {}
   #
   # Read the atom section
   print("Reading atoms section..")
   for ii in range(n_atoms):
      line_split = lines[Atoms_start + ii].split()

      iat = Atom()
      iat.Id    = int(  line_split[DATA_COL_ID])
      iat.molId = int(  line_split[DATA_COL_MOLID])
      iat.type  = int(  line_split[DATA_COL_TYPE])
      iat.q     = float(line_split[DATA_COL_Q])
      iat.x     = float(line_split[DATA_COL_X])
      iat.y     = float(line_split[DATA_COL_Y])
      iat.z     = float(line_split[DATA_COL_Z])
      iat.x     += int(line_split[DATA_COL_IX])*Lx
      iat.y     += int(line_split[DATA_COL_IY])*Ly
      iat.z     += int(line_split[DATA_COL_IZ])*Lz
      atom_list[iat.Id] = deepcopy(iat)
   print()

   print("Coarse graining atomistic trajectories..")

   n_beads = 0
   for group in network:
      for bead in group.beads.values():
         bead_x = 0
         bead_y = 0
         bead_z = 0
         bead_q = 0
         bead_mass = 0
         for aId in bead.aIds:
            atom = atom_list[aId]
            imass = mass_of_type[atom.type]

            bead_q += atom.q
            bead_x += atom.x * imass
            bead_y += atom.y * imass
            bead_z += atom.z * imass
            bead_mass += imass

         bead.x = bead_x / bead_mass
         bead.y = bead_y / bead_mass
         bead.z = bead_z / bead_mass
         bead.q = bead_q
         bead.mass = bead_mass

         n_beads += 1
         #print(bead.bId, bead.molId, num_of_type[bead.type], bead.q, bead.x, bead.y, bead.z," # ",bead.mass, bead.type)

   print("Assigning mass to beads..")
   mass_of_beads = {}
   for group in network:
      for bead in group.beads.values():
         if not bead.type in mass_of_beads:
            mass_of_beads[bead.type] = []
         mass_of_beads[bead.type].append(bead.mass)
   from statistics import mean
   mass_of_bead = {}
   for key in mass_of_beads:
      mass_of_bead[key] = mean(mass_of_beads[key])

   return box, mass_of_type, atom_list, network, mass_of_bead, network

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
      f.write("%s coeff_pair_%s # %s\n" % (num_of_type[bead_type], bead_type, bead_type))

   if num_of_bond_type:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for bond in num_of_bond_type:
         f.write("%s coeff_bond_%s # %s\n" % (num_of_bond_type[bond], bond, bond))

   if num_of_angle_type:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for angle in num_of_angle_type:
         f.write("%s coeff_angle_%s # %s\n" % (num_of_angle_type[angle], angle, angle))

   f.close()



def ExportLammpsDump(path_dump_in, path_dump_out, n_frame, EV_FRAME, network, mass_of_type, num_of_type):
   dump_col = GetDumpFormat(path_dump_in)

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

      if n_frame != "MAX":
         if (iframe > 10 and iframe % int(n_frame / 10.0) == 0):
            print(iframe)
         if iframe == n_frame:
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
         iat = Atom()
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
