import sys

from network import c_bead_bond, c_bead_angle, print_net_stats

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
         print_net_stats(group)

def EnumerateTypes(aLocIds, bond_pairs, angle_pairs):
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

   return num_of_type, num_of_bond_type, num_of_angle_type

def GenBond(network, num_of_bond_type):
   bead_bonds = []
   bead_bondId = 1
   n_beadbonds = 0
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
   return bead_bonds

def GenAngle(network, num_of_angle_type, bead_bonds):
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

   return bead_angles
