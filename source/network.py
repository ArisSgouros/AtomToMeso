import sys

from copy import deepcopy

EMPTY = "empty"
H_MASS = 1.008

class Atom:
   def __init__(self):
      self.Id = EMPTY
      self.molId = EMPTY
      self.type = EMPTY
      self.q = ""
      self.x = ""
      self.y = ""
      self.z = ""

class BeadBond:
   def __init__(self, Id, type, i, j, itype, jtype):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.itype = itype
      self.jtype = jtype

class BeadAngle:
   def __init__(self, Id, type, i, j, k, itype, jtype, ktype):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.k = k
      self.itype = itype
      self.jtype = jtype
      self.ktype = ktype

class Bead():
   def __init__(self, bId, type, aIds):
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

class Group():
   def __init__(self, type):
      self.type = type

      self.beads = {}
      self.nbead = 0
      self.bonds = []
      self.nbond = 0
      self.nat = 0

   def add_bead(self, bId, type, aIds):
      self.beads[bId] = Bead(bId, type, aIds)
      self.nbead += 1
      self.nat += self.beads[bId].nat

   def add_bond(self, ii, jj):
      self.bonds.append( [self.beads[ii], self.beads[jj]] )
      self.nbond += 1

   def add_bead_by_ref(self, bead):
      self.beads[bead.bId] = bead
      self.nbead += 1
      self.nat += bead.nat

   def add_bond_by_ref(self, bond):
      self.bonds.append(bond)
      self.nbond += 1

   def set_bhead(self, bId):
      self.bhead = self.beads[bId]

   def set_btail(self, bId):
      self.btail = self.beads[bId]

   def ins_btail(self, bead):
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

      # generate a deepcopy because Group is imutable
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
         left_group.add_bead_by_ref(bead)
         #left_group.nat += bead.nat

      # Connect the head of the right group to the tail of the left
      left_group.add_bond_by_ref( [left_group.btail, right_group.bhead]  )

      # Append the bonds of right_group to left_group
      for bond in right_group.bonds:
         left_group.add_bond_by_ref(bond)

      left_group.ins_btail(right_group.btail)

      # Generate a new type description
      left_group.type += "-" + right_group.type

      #print(left_group.type)
      #for bead in left_group.beads.values():
      #   print (bead.bId, bead.ahead, bead.atail)

      return left_group


def PrintNetStat(group):
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
