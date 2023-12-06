import sys

sys.path.insert(0, '../../source/')

from copy import deepcopy
from export import ExportLammpsData, ExportLammpsDump
from network import EMPTY, H_MASS, c_atom, c_bead_bond, c_bead_angle, c_bead, c_group, print_net_stats
from process import Process

LAMMPS_DATA_INPUT = "spe8_w850.data"
LAMMPS_DUMP_INPUT = "spe8_w850_0.lammpstrj"
LAMMPS_DATA_OUTPUT = "o.cg.dat"
LAMMPS_DUMP_OUTPUT = "o.cg.lammpstrj"
verbose = True
N_FRAME = 40
EV_FRAME = 1
atomtype = 'full'

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

# Print info regarding each group
groups = [gSPE, gW]

print("-----------------------------------------------")
print("Base group stats")
print("-----------------------------------------------")
for group in groups:
   print_net_stats(group)
print()

# Define the bond types
bond_pairs = [ ["SB","SB"], ["SO","SB"], ["SO","SN"], ["SN","SS"] ]

# Define the angle types
angle_pairs = [ ["SB","SB","SB"], ["SB","SB","SO"], ["SB","SO","SN"], ["SO","SN","SS"] ]

# Report IDs assigned to beads
print("-----------------------------------------------")
print("Local IDs of beads")
print("-----------------------------------------------")
print("type len LocalIds")
for key in aLocIds:
   print(key, len(aLocIds[key]), aLocIds[key])
print()

network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box, mass_of_type = Process(aLocIds, bond_pairs, angle_pairs, network, LAMMPS_DATA_INPUT, atomtype, verbose)

print("-----------------------------------------------")
print("Generating the new lammps datafile with name " + LAMMPS_DATA_OUTPUT + "..")
print("-----------------------------------------------")
ExportLammpsData(LAMMPS_DATA_OUTPUT, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box)
print("SUCCESS!")

print("-----------------------------------------------")
print("Generating the new lammps dump with name " + LAMMPS_DUMP_OUTPUT + "..")
print("-----------------------------------------------")
ExportLammpsDump(LAMMPS_DUMP_INPUT, LAMMPS_DUMP_OUTPUT, N_FRAME, EV_FRAME, network, mass_of_type, num_of_type)
