##!/bin/bash

# run msi2lmp
/media/aps/D/Programs/msi2lmp_aps/src/msi2lmp.exe dgebf_a30 -c 0 -f oplsaa -print 3 -nocenter -dir_frc source/. -skipcoeff #> o.log

# convert to lammpstrj
python /media/aps/D/Programs/LmpTool/DataToDump/data_to_dump.py dgebf_a30.data full -file_type="lammpstrj" --fmt="id,mol,type,x,y,z" -dump_file="o.dump.lammpstrj"

# coarse-grain
python atom_to_meso.py #> o.log

# generate the SMART file
#python source/lmp_gen_SMART_MOL_v5.py smart_files 1 DGEBF

# apply the OPLS force field
#python source/lmp_apply_OPLS_v4.py dgebf.data SMART_1-DGEBF.dat source/oplsaa_from_foyer_group_ORIG.xml
