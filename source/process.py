import sys

from network import c_bead_bond

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

