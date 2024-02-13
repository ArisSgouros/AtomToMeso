import sys

# import path of source files from command line
try:
   path_a2m=sys.argv[1]
except:
   path_a2m='../../source/'
sys.path.insert(0, path_a2m)

from copy import deepcopy
from export import ExportLammpsData, ExportLammpsDump
from network import Atom, BeadBond, BeadAngle, Bead, Group, PrintNetStat, ReadGroup, ExportGroup
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

   gSPE = ReadGroup("SPE")
 
   # text the export function
   gSPE_COPY = deepcopy(gSPE)
   gSPE_COPY.type = "SPE_COPY"
   ExportGroup(gSPE_COPY)

   gW = ReadGroup("W")
   gW_COPY = deepcopy(gW)
   gW_COPY.type = "W_COPY"
   ExportGroup(gW_COPY)

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

   print("-----------------------------------------------")
   print("Testing the production section..")
   print("-----------------------------------------------")
   print('mean, std of bond 1  : ', ComputeBond(path_lammps_data_out, [1], n_frame, path_lammps_dump_out))
   print('mean, std of bond 2  : ', ComputeBond(path_lammps_data_out, [2], n_frame, path_lammps_dump_out))
   print('mean, std of bond 3  : ', ComputeBond(path_lammps_data_out, [3], n_frame, path_lammps_dump_out))
   print('mean, std of bond 4  : ', ComputeBond(path_lammps_data_out, [4], n_frame, path_lammps_dump_out))
   print('mean, std of angle 1 : ', ComputeAngle(path_lammps_data_out, [1], n_frame, path_lammps_dump_out))
   print('mean, std of angle 2 : ', ComputeAngle(path_lammps_data_out, [2], n_frame, path_lammps_dump_out))
   print('mean, std of angle 3 : ', ComputeAngle(path_lammps_data_out, [3], n_frame, path_lammps_dump_out))
   print('mean, std of angle 4 : ', ComputeAngle(path_lammps_data_out, [4], n_frame, path_lammps_dump_out))
