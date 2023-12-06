import sys

from copy import deepcopy


def ExportLammpsData(filename, network, bead_bonds, bead_angles, mass_of_bead, num_of_bond_type, num_of_angle_type):
   f = open(filename, 'w')

   f.write('# Coarse grained LAMMPS data file\n')
   f.write('\n')
   f.write('%d atoms\n' % n_beads)
   f.write('%d bonds\n' % n_beadbonds)
   f.write('%d angless\n' % n_beadangles)
   f.write('\n')
   f.write('%d atom types\n' % n_bead_types)
   f.write('%d bond types\n' % n_bead_bondtypes)
   f.write('%d angle types\n' % n_bead_angletypes)
   f.write('\n')
   f.write('%s %s xlo xhi\n' % (xlo, xhi))
   f.write('%s %s ylo yhi\n' % (ylo, yhi))
   f.write('%s %s zlo zhi\n' % (zlo, zhi))
   f.write('\n')
   f.write('Atoms\n')
   f.write('\n')

   for group in network:
      for bead in group.beads.values():
         f.write("%d %d %d %16.9f %16.9f %16.9f %16.9f # %s\n" % ( bead.bId, bead.molId, num_of_type[bead.type], bead.q, bead.x, bead.y, bead.z, bead.type ))

   if n_beadbonds:
      f.write('\n')
      f.write('Bonds\n')
      f.write('\n')
      for bond in bead_bonds:
         f.write("%d %d %d %d # %s_%s\n" % ( bond.Id, bond.type, bond.i, bond.j, bond.itype, bond.jtype ))

   if n_beadangles:
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

   if n_beadbonds:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for bond in num_of_bond_type:
         f.write(str(num_of_bond_type[bond])+" Bond coeffs # " + bond+"\n")

   if n_beadangles:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for angle in num_of_angle_type:
         f.write(str(num_of_angle_type[angle])+" Angle coeffs # " + angle+"\n")

   f.close()



def ExportLammpsDump(path_dump_in, path_dump_out, N_FRAME, EV_FRAME, network):
   fin = open(path_dump_in, 'r')
   fout = open(path_dump_out, 'w')

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




n_sysargv = len(sys.argv)
if n_sysargv != 2 and n_sysargv != 5:
   print
   print ("The required formats are the following:")
   print
   print ("python atom_to_meso.py \"LAMMPS_DATA_INPUT\" ")
   print
   print ("example: python atom_to_meso.py SVSV-W800.dat")
   print
   print ("or")
   print
   print ("python atom_to_meso.py \"LAMMPS_DATA_INPUT\" \"LAMMPS_DUMP_INPUT\" \"N_FRAME\" \"EV_FRAME\"")
   print ("example 1: python atom_to_meso.py SS-W800.dat 1_prod.lammpstrj MAX 5")
   print ("example 2: python atom_to_meso.py SS-W800.dat 1_prod.lammpstrj 30 10")
   print
   print ("exiting..")
   sys.exit()

LAMMPS_DATA_INPUT = sys.argv[1]
LAMMPS_DATA_OUTPUT = "o.cg.dat"

if n_sysargv == 5:
   LAMMPS_DUMP_INPUT = sys.argv[2]

   N_FRAME = sys.argv[3]
   if N_FRAME != "MAX":
      N_FRAME = int(sys.argv[3])

   EV_FRAME = int(sys.argv[4])

   LAMMPS_DUMP_OUTPUT = "o.cg.lammpstrj"

EMPTY = "empty"
H_MASS = 1.008

class c_atom:
   def __init__(self):
      self.Id = EMPTY
      self.molId = EMPTY
      self.type = EMPTY
      self.q = ""
      self.x = ""
      self.y = ""
      self.z = ""

class c_bead_bond:
   def __init__(self, Id, type, i, j, itype, jtype):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.itype = itype
      self.jtype = jtype

class c_bead_angle:
   def __init__(self, Id, type, i, j, k, itype, jtype, ktype):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.k = k
      self.itype = itype
      self.jtype = jtype
      self.ktype = ktype

class c_bead():
   def __init__(self, type, aIds, bId):
      self.type = type
      self.molId = EMPTY
      self.bId = bId
      self.aIds = aIds
      self.nat = len(aIds)
      self.ahead = EMPTY
      self.atail = EMPTY
      self.x = ""
      self.y = ""
      self.z = ""
      self.q = ""

   def set_ahead(self, aId):
      self.ahead = aId

   def set_atail(self, aId):
      self.atail = aId

   def shift_aIds(self, shift):
      self.aIds = [ aId + shift for aId in self.aIds ]
      #print(bead.bId, bead.aIds, shift)

      if self.ahead != EMPTY:
         self.ahead += shift
      if self.atail != EMPTY:
         self.atail += shift

class c_group():
   def __init__(self, type):
      self.type = type

      self.beads = {}
      self.nbead = 0
      self.bonds = []
      self.nbond = 0
      self.nat = 0

   def add_bead(self, bead):
      self.beads[bead.bId] = bead
      self.nbead += 1
      self.nat += bead.nat

   def add_bond(self, bond):
      self.bonds.append(bond)
      self.nbond += 1

   def set_bhead(self, bead):
      self.bhead = bead

   def set_btail(self, bead):
      self.btail = bead

   def rmv_atom_from_bead(self, bead, aId):
      bead.aIds.remove(aId)
      bead.nat -= 1
      self.nat -= 1

   def shift_lids(self, shift):
      aux = {}
      for bead in self.beads.values():
         bead.bId += shift
         aux[bead.bId] = bead
      self.beads = deepcopy(aux)

   # merge groups
   def __add__(self, new_group):

      # generate a deepcopy because c_group is imutable
      right_group = deepcopy(new_group)
      left_group = deepcopy(self)

      # increment the bead ids of the right group
      left_group_nbead = left_group.nbead
      for bead in right_group.beads.values():
         bead.bId += left_group.nbead

      # pop head atoms from the head beads of the right group
      for bead in right_group.beads.values():
        if bead is right_group.bhead:
           #print('AHEAD BF: ahead ',bead.ahead," from ", bead.aIds, "( bead nat: ", bead.nat , " gr nat: ", right_group.nat)
           right_group.rmv_atom_from_bead(bead, bead.ahead)
           bead.ahead = EMPTY
           #print('AHEAD AF: ahead ',bead.ahead," from ", bead.aIds, "( bead nat: ", bead.nat , " gr nat: ", right_group.nat)

      # pop tail atoms from the tail beads of the left group
      for bead in left_group.beads.values():
        if bead is left_group.btail:
           #print('ATAIL BF: atail ',bead.atail," from ", bead.aIds, "( bead nat: ", bead.nat , " gr nat: ", left_group.nat)
           left_group.rmv_atom_from_bead(bead, bead.atail)
           bead.atail = EMPTY
           #print('ATAIL AF: atail ',bead.atail," from ", bead.aIds, "( bead nat: ", bead.nat , " gr nat: ", left_group.nat)
      #print(left_group.type," + ",right_group.type)
      #print("left group n atoms: "+str(left_group.nat) + " " + str(right_group.nat))

      shift = left_group.nat - 1
      #print("SHIFT: "+str(shift))

      # Update and increment the ids of the right group atoms
      for bead in right_group.beads.values():
         #print(bead.bId, bead.aIds)
         bead.shift_aIds(shift)

      # Append the beads of right_group to left_group
      for bead in right_group.beads.values():
         left_group.add_bead(bead)
         #left_group.nat += bead.nat

      # Connect the head of the right group to the tail of the left
      left_group.add_bond( [left_group.btail, right_group.bhead]  )

      # Append the bonds of right_group to left_group
      for bond in right_group.bonds:
         left_group.add_bond(bond)

      left_group.set_btail(right_group.btail)

      # Generate a new type description
      left_group.type += "-" + right_group.type

      #print(left_group.type)
      #for bead in left_group.beads.values():
      #   print (bead.bId, bead.ahead, bead.atail)

      return left_group


def print_net_stats(group):
   print('---',group.type,'---')
   print("nbeads: ", group.nbead)
   print("bIds type ahead atail aIds")
   try:
      for bead in group.beads.values():
         print(bead.bId, bead.type, bead.ahead, bead.atail, bead.aIds)
      print("head",group.bhead.bId, group.bhead.type)
      print("tail",group.btail.bId, group.btail.type)
   except:
      print("NO HEADS/TAILS")
   if group.bonds:
      print("bond pairs")
      for bond in group.bonds:
         print(bond[0].bId, bond[1].bId, " # ", bond[0].type, bond[1].type)
   else:
      print("NO BONDS")
   print()

############################################################
# Interface
############################################################

# Set the format of Lammps data and dump files
DUMP_COL_ID    = 0
DUMP_COL_MOLID = 1
DUMP_COL_TYPE  = 2
DUMP_COL_X     = 3
DUMP_COL_Y     = 4
DUMP_COL_Z     = 5

DATA_COL_ID    = 0
DATA_COL_MOLID = 1
DATA_COL_TYPE  = 2
DATA_COL_Q     = 3
DATA_COL_X     = 4
DATA_COL_Y     = 5
DATA_COL_Z     = 6
DATA_COL_IX    = DATA_COL_X + 3
DATA_COL_IY    = DATA_COL_Y + 3
DATA_COL_IZ    = DATA_COL_Z + 3

# Set the local ids for each segment represented by a bead
aLocIds={}
aLocIds["SB"]=[1,2,3,4,5,6,38,39,40,41]
aLocIds["SO"]=[7,8,9,10,11,12]
aLocIds["SN"]=[13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
aLocIds["SS"]=[28,29,30,31,32,33,34,35,36,37]
aLocIds["W"]=[1,2,3]

# Create the bead objects
# SPE beads
bSB = c_bead("SB",aLocIds["SB"],1)
bSB.set_ahead(1)
bSB.set_atail(41)
bSO = c_bead("SO",aLocIds["SO"],2)
bSN = c_bead("SN",aLocIds["SN"],3)
bSS = c_bead("SS",aLocIds["SS"],4)

# WATER bead
bW = c_bead("W",aLocIds["W"],1)

# Create SPE group (repeat unit)
gSPE = c_group("SPE")
gSPE.add_bead( bSB )
gSPE.add_bead( bSO )
gSPE.add_bead( bSN )
gSPE.add_bead( bSS )
gSPE.add_bond( [bSB, bSO] )
gSPE.add_bond( [bSO, bSN] )
gSPE.add_bond( [bSN, bSS] )
gSPE.set_bhead(bSB)
gSPE.set_btail(bSB)

# Create the water group
gW = c_group("W")
gW.add_bead(bW)

# Generate chains by merging groups
print("Generating chains..")
S8 = gSPE+gSPE+gSPE+gSPE+gSPE+gSPE+gSPE+gSPE

# generate a network with 850 water and 1 SPE8 molecules
network = []
for ii in range(850):
   network.append(deepcopy(gW))
network.append(deepcopy(S8))

# Define the bond types
bond_pairs = [ ["SB","SB"], ["SO","SB"], ["SO","SN"], ["SN","SS"] ]

# Define the angle types
angle_pairs = [ ["SB","SB","SB"], ["SB","SB","SO"], ["SB","SO","SN"], ["SO","SN","SS"] ]

########################################
# Internal Section
########################################

# Report IDs assigned to beads
print("-----------------------------------------------")
print("Local IDs of beads")
print("-----------------------------------------------")
print("type len LocalIds")
for key in aLocIds:
   print(key, len(aLocIds[key]), aLocIds[key])
print()

# Enumerate bead types
id = 1
num_of_type = {}
for key in aLocIds:
   num_of_type[key] = id
   id += 1

# get the number of bead types
aux = set()
for type in num_of_type.values():
   aux.add(type)
n_bead_types = len(aux)
print(n_bead_types)

# Enumerate the bond types
num_of_bond_type = {}
n_bead_bondtypes = len(bond_pairs)
id = 1
for bond_pair in bond_pairs:
   if bond_pair[0] > bond_pair[-1]:
      bond_pair.reverse()
   aux = "_".join(bond_pair)
   num_of_bond_type[aux] = id
   id += 1

# Enumerate the angle types
num_of_angle_type = {}
n_bead_angletypes = len(angle_pairs)
id = 1
for angle_pair in angle_pairs:
   if angle_pair[0] > angle_pair[-1]:
      angle_pair.reverse()
   aux = "_".join(angle_pair)
   num_of_angle_type[aux] = id
   id += 1

# Print info regarding each group
groups = [gSPE, gW]

print("-----------------------------------------------")
print("Base group stats")
print("-----------------------------------------------")
for group in groups:
   print_net_stats(group)
print()


print("-----------------------------------------------")
print("Generating global IDs & Network stats..")
print("-----------------------------------------------")
gbId = 0
gaId = 0
molId = 1
for group in network:
   group.shift_lids(gbId)
   for bead in group.beads.values():
      bead.shift_aIds(gaId)
      bead.molId = molId
   gbId += group.nbead
   gaId += group.nat
   molId += 1
   print_net_stats(group)

bead_bonds = []
bead_bondId = 1
n_beadbonds = 0
print("-----------------------------------------------")
print("Generating bonds between beads")
print("-----------------------------------------------")
for group in network:
   for bond in group.bonds:

      aux = "_".join(sorted([bond[0].type,bond[1].type]))
      try:
         type = num_of_bond_type[ aux ]
      except:
         print("Unidentified bond type: " + aux)
         exit()

      i = bond[0].bId
      j = bond[1].bId
      itype = bond[0].type
      jtype = bond[1].type

      ibond = c_bead_bond(bead_bondId,type,i,j,itype,jtype)
      bead_bonds.append(ibond)
      bead_bondId += 1
      n_beadbonds += 1


print("-----------------------------------------------")
print("Generating angles between subsequent beads")
print("-----------------------------------------------")
bead_angles = []
n_beadangles = 0

all_angle_ids = []

angle_id = 1
for abond in bead_bonds:
   for bbond in bead_bonds:
      if abond is bbond: continue

      angle_ids = []
      if abond.i is bbond.i:
         angle_ids = [abond.j, abond.i, bbond.j]
         angle_types = [abond.jtype, abond.itype, bbond.jtype]
      if abond.i is bbond.j:
         angle_ids = [abond.j, abond.i, bbond.i]
         angle_types = [abond.jtype, abond.itype, bbond.itype]
      if abond.j is bbond.i:
         angle_ids = [abond.i, abond.j, bbond.j]
         angle_types = [abond.itype, abond.jtype, bbond.jtype]
      if abond.j is bbond.j:
         angle_ids = [abond.i, abond.j, bbond.i]
         angle_types = [abond.itype, abond.jtype, bbond.itype]

      if angle_ids:
         # Check if the angle already exists
         aux_sort_ids = sorted(angle_ids)
         if aux_sort_ids in all_angle_ids:
            continue
         else:
            all_angle_ids.append(aux_sort_ids)

         # Sort the first and last types for them to be unique
         aux_type = [ii for ii in angle_types]
         if aux_type[0] > aux_type[-1]:
            aux_type.reverse()

         aux_type = "_".join(aux_type)

         try:
            type = num_of_angle_type[ aux_type ]
         except:
            print("Unidentified angle type: " + aux_type)
          #  exit()

         iangle = c_bead_angle(angle_id, type, angle_ids[0], angle_ids[1], angle_ids[2], angle_types[0], angle_types[1], angle_types[2])
         bead_angles.append(iangle)


         angle_id += 1
         n_beadangles += 1

print("-----------------------------------------------")
print("Generating coarse grained LAMMPS DATA FILE"     )
print("-----------------------------------------------")
print()
print("Parsing lammps data file: ",LAMMPS_DATA_INPUT)
print()
n_atoms = 0
g = open(LAMMPS_DATA_INPUT, 'r')

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
      xlo = line_split[0]
      xhi = line_split[1]
      Lx = float(xhi)-float(xlo)
   if "ylo yhi" in line:
      ylo = line_split[0]
      yhi = line_split[1]
      Ly = float(yhi)-float(ylo)
   if "zlo zhi" in line:
      zlo = line_split[0]
      zhi = line_split[1]
      Lz = float(zhi)-float(zlo)

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

   iat = c_atom()
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

mass_of_bead = {}

print("Assigning mass to beads..")
for group in network:
   for bead in group.beads.values():
      imass = deepcopy(bead.mass)
      if not bead.type in mass_of_bead:

         print(bead.type)
         print(imass)
         if bead.ahead is not EMPTY:
            imass -= H_MASS
         print(imass, bead.ahead)
         if bead.atail is not EMPTY:
            imass -= H_MASS
         print(imass, bead.atail)

         mass_of_bead[bead.type] = imass





print("-----------------------------------------------")
print("Generating the new lammps datafile with name " + LAMMPS_DATA_OUTPUT + "..")
print("-----------------------------------------------")
ExportLammpsData(LAMMPS_DATA_OUTPUT, network, bead_bonds, bead_angles, mass_of_bead, num_of_bond_type, num_of_angle_type)
print("SUCCESS!")



if n_sysargv == 5:
   print("-----------------------------------------------")
   print("Generating the new lammps dump with name " + LAMMPS_DUMP_OUTPUT + "..")
   print("-----------------------------------------------")
   ExportLammpsDump(LAMMPS_DUMP_INPUT, LAMMPS_DUMP_OUTPUT, N_FRAME, EV_FRAME, network)
