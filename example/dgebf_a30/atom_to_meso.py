import sys

sys.path.insert(0, '../../source/')

from copy import deepcopy
from export import ExportLammpsData, ExportLammpsDump
from network import Atom, BeadBond, BeadAngle, Bead, Group, PrintNetStat
from process import Process

if __name__ == "__main__":
   path_lammps_data_in = "dgebf_a30_mod.data"
   path_lammps_dump_in = "o.dump.lammpstrj"
   path_lammps_data_out = "o.cg.dat"
   path_lammps_dump_out = "o.cg.lammpstrj"
   verbose = True
   n_frame = 1
   ev_frame = 1
   atomtype = 'full'

   gDGEBF = Group("DGFBE")
   gDGEBF.add_bead(1,"GR",[1,44,5,2,7,4,3,47,6,9,10])
   gDGEBF.add_bead(2,"GO",[8])
   gDGEBF.add_bead(3,"GP",[11,12,13,17,15,18,14,16,19,20])
   gDGEBF.add_bead(4,"GM",[21,42,43])
   gDGEBF.add_bead(5,"GP",[22,23,24,28,26,27,29,25,31,30])
   gDGEBF.add_bead(6,"GO",[32])
   gDGEBF.add_bead(7,"GR",[36,35,33,46,37,34,39,40,38,45,41])

   gDGEBF.add_bond(1,2)
   gDGEBF.add_bond(2,3)
   gDGEBF.add_bond(3,4)
   gDGEBF.add_bond(4,5)
   gDGEBF.add_bond(5,6)
   gDGEBF.add_bond(6,7)

   gDGEBF.set_bhead(7)
   gDGEBF.set_btail(1)
   gDGEBF.beads[7].set_ahead(45)
   gDGEBF.beads[1].set_atail(1)



   # Generate chains by merging groups
   print("Generating chains..")

   # generate a network with 850 water and 1 SPE8 molecules
   network = [deepcopy(gDGEBF) for i in range(30)]

   # Print info regarding each group
   groups = [gDGEBF]

   print("-----------------------------------------------")
   print("Base group stats")
   print("-----------------------------------------------")
   for group in groups:
      PrintNetStat(group)
   print()

   # Report IDs assigned to beads
   print("-----------------------------------------------")
   print("Local IDs of beads")
   print("-----------------------------------------------")
   print("type len LocalIds")
   for bead in gDGEBF.beads.values():
      print(bead.type, len(bead.aIds), bead.aIds)
   print()

   network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box, mass_of_type = Process(network, path_lammps_data_in, atomtype, verbose)

   print("-----------------------------------------------")
   print("Generating the new lammps datafile with name " + path_lammps_data_out + "..")
   print("-----------------------------------------------")
   ExportLammpsData(path_lammps_data_out, network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box)

   print("-----------------------------------------------")
   print("Generating the new lammps dump with name " + path_lammps_dump_out + "..")
   print("-----------------------------------------------")
   ExportLammpsDump(path_lammps_dump_in, path_lammps_dump_out, n_frame, ev_frame, network, mass_of_type, num_of_type)
