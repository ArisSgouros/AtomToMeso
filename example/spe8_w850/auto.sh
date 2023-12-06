DATAFILE="spe8_w850.data"
DUMP_FILE="spe8_w850_0.lammpstrj"
n_frame=40
EVFRAME=1

python atom_to_meso.py > o.log

#python ../../source/lmp_type_strip.py o.cg.lammpstrj 3 "[5]" > o.cg.strip_5.lammpstrj
#
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [1] $n_frame o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [2] $n_frame o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [3] $n_frame o.cg.lammpstrj
#python ../../source/lmp_angle_dist_v3.py o.cg.dat [4] $n_frame o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [1] $n_frame o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [2] $n_frame o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [3] $n_frame o.cg.lammpstrj
#python ../../source/lmp_bond_dist_v3.py o.cg.dat [4] $n_frame o.cg.lammpstrj
