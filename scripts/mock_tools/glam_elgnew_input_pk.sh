#!/bin/bash

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 4 clustering-stats --tracer ELG_input --zrange 0.8 1.1 1.1 1.6 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/glam_new/mock$i/ --stats_dir /pscratch/sd/a/ajross/glam_new/mock$i/ --combine --nran 18 --boxsize 10000 --cellsize 10

done