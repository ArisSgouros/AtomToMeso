import sys

sys.path.insert(0, '../../source/')

from copy import deepcopy
from export import ExportLammpsData, ExportLammpsDump
from network import EMPTY, H_MASS, Atom, BeadBond, BeadAngle, Bead, Group, PrintNetStat
from process import Process

if __name__ == "__main__":
   path_lammps_data_in = "spe8_w850.data"
   path_lammps_dump_in = "spe8_w850_0.lammpstrj"
   path_lammps_data_out = "o.cg.dat"
   path_lammps_dump_out = "o.cg.lammpstrj"
   verbose = True
   n_frame = 40
   ev_frame = 1
   atomtype = 'full'

   # Create the bead objects
   # SPE beads
   bSB = Bead(1,"SB",[1,2,3,4,5,6,38,39,40,41])
   bSB.set_ahead(1)
   bSB.set_atail(41)
   bSO = Bead(2,"SO",[7,8,9,10,11,12])
   bSN = Bead(3,"SN",[13,14,15,16,17,18,19,20,21,22,23,24,25,26,27])
   bSS = Bead(4,"SS",[28,29,30,31,32,33,34,35,36,37])

   # WATER bead
   bW = Bead(1,"W",[1,2,3])

   bead_types = [bSB, bSO, bSN, bSS, bW]

   # Create SPE group (repeat unit)
   gSPE = Group("SPE")
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
   gW = Group("W")
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
      PrintNetStat(group)
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
   for bead in bead_types:
      print(bead.type, len(bead.aIds), bead.aIds)
   print()

   network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box, mass_of_type = Process(bead_types, bond_pairs, angle_pairs, network, path_lammps_data_in, atomtype, verbose)

   print("-----------------------------------------------")
   print("Generating the new lammps datafile with name " + path_lammps_data_out + "..")
   print("-----------------------------------------------")
   ExportLammpsData(path_lammps_data_out, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box)
   print("SUCCESS!")

   print("-----------------------------------------------")
   print("Generating the new lammps dump with name " + path_lammps_dump_out + "..")
   print("-----------------------------------------------")
   ExportLammpsDump(path_lammps_dump_in, path_lammps_dump_out, n_frame, ev_frame, network, mass_of_type, num_of_type)
