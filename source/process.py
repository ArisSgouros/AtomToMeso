import sys

from network import BeadBond, BeadAngle, PrintNetStat
from export import ParseLammpsData

def GenGlobalId(network, verbose):
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
      if verbose:
         PrintNetStat(group)

def EnumerateTypes(network):
   itype_max = 0
   num_of_type = {}
   for group in network:
      for bead in group.beads.values():
         aux = bead.type
         if not aux in num_of_type:
            itype_max += 1
            num_of_type[aux] = itype_max
   return num_of_type

def GenBond(network):
   bead_bonds = []
   bondId = 1
   btype_max = 0
   num_of_bond_type = {}
   for group in network:
      for bond in group.bonds:

         aux = "_".join(sorted([bond[0].type,bond[1].type]))
         if not aux in num_of_bond_type:
            btype_max += 1
            num_of_bond_type[aux] = btype_max
         type = num_of_bond_type[ aux ]

         i = bond[0].bId
         j = bond[1].bId
         itype = bond[0].type
         jtype = bond[1].type

         ibond = BeadBond(bondId,type,i,j,itype,jtype)
         bead_bonds.append(ibond)
         bondId += 1
   return bead_bonds, num_of_bond_type

def GenAngle(network, bead_bonds):
   bead_angles = []

   all_angle_ids = []

   atype_max = 0
   num_of_angle_type = {}

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
            if not aux_type in num_of_angle_type:
               atype_max += 1
               num_of_angle_type[aux_type] = atype_max
            type = num_of_angle_type[ aux_type ]

            iangle = BeadAngle(angle_id, type, angle_ids[0], angle_ids[1], angle_ids[2], angle_types[0], angle_types[1], angle_types[2])
            bead_angles.append(iangle)


            angle_id += 1

   return bead_angles, num_of_angle_type

def Process(network, path_lammps_data_in, atom_type, verbose):

   num_of_type = EnumerateTypes(network)

   print("-----------------------------------------------")
   print("Generating global IDs & Network stats..")
   print("-----------------------------------------------")
   GenGlobalId(network, verbose)

   # Check if the atom ids are unique
   aux_id_list = []
   for group in network:
      for bead in group.beads.values():
         for aid in bead.aIds:
            print(aid)
            if aid in aux_id_list:
               print("ERROR: nonunique atom id", aid)
               sys.exit()
            else:
               aux_id_list.append(aid)

   print("-----------------------------------------------")
   print("Generating bonds between beads")
   print("-----------------------------------------------")
   bead_bonds, num_of_bond_type = GenBond(network)

   print("-----------------------------------------------")
   print("Generating angles between subsequent beads")
   print("-----------------------------------------------")
   bead_angles, num_of_angle_type = GenAngle(network, bead_bonds)

   print("-----------------------------------------------")
   print("Generating coarse grained LAMMPS DATA FILE"     )
   print("-----------------------------------------------")
   print()
   print("Parsing lammps data file: ",path_lammps_data_in)
   print()

   box, mass_of_type, atom_list, network, mass_of_bead, network = ParseLammpsData(path_lammps_data_in, network, atom_type)

   return network, bead_bonds, bead_angles, mass_of_bead, num_of_type, num_of_bond_type, num_of_angle_type, box, mass_of_type



