import sys

# import path of source files from command line
try:
   path_a2m=sys.argv[1]
except:
   path_a2m='../../source/'
sys.path.insert(0, path_a2m)

from copy import deepcopy
from export import ExportLammpsData, ExportLammpsDump
from network import Atom, BeadBond, BeadAngle, Bead, Group, PrintNetStat
from process import Process
from compute_bond import ComputeBond
from compute_angle import ComputeAngle

if __name__ == "__main__":
   path_lammps_data_in = "spe8_w850.data"
   path_lammps_dump_in = "spe8_w850_0.lammpstrj"
   path_lammps_data_out = "o.cg.dat"
   path_lammps_dump_out = "o.cg.lammpstrj"
   verbose = True
   n_frame = 40
   ev_frame = 1
   atomtype = 'full'

   # Create SPE group (repeat unit)
   gSPE = Group("SPE")
   gSPE.add_bead( 1, "SB", [1,2,3,4,5,6,38,39,40,41] )
   gSPE.add_bead( 2, "SO", [7,8,9,10,11,12] )
   gSPE.add_bead( 3, "SN", [13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] )
   gSPE.add_bead( 4, "SS", [28,29,30,31,32,33,34,35,36,37] )
   gSPE.add_bond( 1, 2 )
   gSPE.add_bond( 2, 3 )
   gSPE.add_bond( 3, 4 )
   gSPE.beads[1].set_ahead(41)
   gSPE.beads[1].set_atail(1)
   gSPE.set_bhead(1)
   gSPE.set_btail(1)

   # Create the water group
   gW = Group("W")
   gW.add_bead(1,"W",[1,2,3])

   # Generate chains by merging groups
   print("Generating chains..")
   S8 = gSPE+gSPE+gSPE+gSPE+gSPE+gSPE+gSPE+gSPE

   # generate a network with 850 water and 1 SPE8 molecules
   network = []
   for ii in range(850):
      network.append(deepcopy(gW))
   network.append(deepcopy(S8))

   # Print info regarding each group
   print("-----------------------------------------------")
   print("Base group stats")
   print("-----------------------------------------------")
   for group in [gSPE, gW]:
      PrintNetStat(group)
   print()

   # Report IDs assigned to beads
   print("-----------------------------------------------")
   print("Local IDs of beads")
   print("-----------------------------------------------")
   print("type len LocalIds")
   for bead in [gSPE.beads[1], gSPE.beads[2], gSPE.beads[3], gSPE.beads[4], gW.beads[1]]:
      print(bead.type, len(bead.aIds), bead.aIds)
   print()

   network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box, mass_of_type = Process(network, path_lammps_data_in, atomtype, verbose)

   print("-----------------------------------------------")
   print("Generating the new lammps datafile with name " + path_lammps_data_out + "..")
   print("-----------------------------------------------")
   ExportLammpsData(path_lammps_data_out, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box)
   print("SUCCESS!")

   print("-----------------------------------------------")
   print("Generating the new lammps dump with name " + path_lammps_dump_out + "..")
   print("-----------------------------------------------")
   ExportLammpsDump(path_lammps_dump_in, path_lammps_dump_out, n_frame, ev_frame, network, mass_of_type, num_of_type)

   ComputeBond(path_lammps_data_out, [1], n_frame, path_lammps_dump_out, lbin=0.006525, verbose=True, debug=True)
   ComputeBond(path_lammps_data_out, [2], n_frame, path_lammps_dump_out, lbin=0.006525, verbose=True, debug=True)
   ComputeBond(path_lammps_data_out, [3], n_frame, path_lammps_dump_out, lbin=0.006525, verbose=True, debug=True)
   ComputeBond(path_lammps_data_out, [4], n_frame, path_lammps_dump_out, lbin=0.006525, verbose=True, debug=True)
   ComputeAngle(path_lammps_data_out, [1], n_frame, path_lammps_dump_out, lbin=0.01, verbose=True, debug=True)
   ComputeAngle(path_lammps_data_out, [2], n_frame, path_lammps_dump_out, lbin=0.01, verbose=True, debug=True)
   ComputeAngle(path_lammps_data_out, [3], n_frame, path_lammps_dump_out, lbin=0.01, verbose=True, debug=True)
   ComputeAngle(path_lammps_data_out, [4], n_frame, path_lammps_dump_out, lbin=0.01, verbose=True, debug=True)
