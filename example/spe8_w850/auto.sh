DATAFILE="spe8_w850.data"
DUMP_FILE="spe8_w850_0.lammpstrj"
path_a2m='../../source'
n_frame=40
EVFRAME=1

python atom_to_meso.py $path_a2m > o.log

#python $path_a2m/compute_bond.py o.cg.dat [1] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_bond.py o.cg.dat [2] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_bond.py o.cg.dat [3] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_bond.py o.cg.dat [4] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_angle.py o.cg.dat [1] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_angle.py o.cg.dat [2] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_angle.py o.cg.dat [3] $n_frame o.cg.lammpstrj
#python $path_a2m/dist_angle.py o.cg.dat [4] $n_frame o.cg.lammpstrj
