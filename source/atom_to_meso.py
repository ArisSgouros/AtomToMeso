import sys

from copy import deepcopy
from export import ParseLammpsData, ExportLammpsData, ExportLammpsDump
from network import EMPTY, H_MASS, c_atom, c_bead_bond, c_bead_angle, c_bead, c_group, print_net_stats
from process import EnumerateTypes, GenBond


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

# Set the format of Lammps data and dump files
dump_col = {'id': 0, 'molid': 1, 'type': 2, 'x': 3, 'y': 4, 'z': 5}
data_col = {'id': 0, 'molid': 1, 'type': 2, 'q': 3, 'x':4, 'y':5, 'z':6, 'ix':7, 'iy':8, 'iz':9}


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

num_of_type, num_of_bond_type, num_of_angle_type = EnumerateTypes(aLocIds, bond_pairs, angle_pairs)



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


print("-----------------------------------------------")
print("Generating bonds between beads")
print("-----------------------------------------------")
bead_bonds = GenBond(network, num_of_bond_type)

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

box, mass_of_type, atom_list, network, mass_of_bead, network = ParseLammpsData(LAMMPS_DATA_INPUT, network, data_col)



print("-----------------------------------------------")
print("Generating the new lammps datafile with name " + LAMMPS_DATA_OUTPUT + "..")
print("-----------------------------------------------")
ExportLammpsData(LAMMPS_DATA_OUTPUT, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box)
print("SUCCESS!")



if n_sysargv == 5:
   print("-----------------------------------------------")
   print("Generating the new lammps dump with name " + LAMMPS_DUMP_OUTPUT + "..")
   print("-----------------------------------------------")
   ExportLammpsDump(LAMMPS_DUMP_INPUT, LAMMPS_DUMP_OUTPUT, N_FRAME, EV_FRAME, network, mass_of_type, dump_col, num_of_type)
