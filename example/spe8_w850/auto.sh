DATAFILE="spe8_w850.data"
DUMP_FILE="spe8_w850_0.lammpstrj"
NFRAME=40
EVFRAME=1

python ../../source/atom_to_meso.py $DATAFILE $DUMP_FILE $NFRAME $EVFRAME > o.log

#python ../../source/lmp_type_strip.py o.cg.lammpstrj 3 "[5]" > o.cg.strip_5.lammpstrj
#
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [1] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [2] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [3] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [4] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [1] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [2] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [3] $NFRAME o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [4] $NFRAME o.cg.lammpstrj
