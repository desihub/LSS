#!/bin/bash

source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 4 clustering-stats --tracer LRG_input --zrange 0.665 0.74 0.74 0.82 0.9 0.9 0.98 0.98 1.1 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --stats_dir /pscratch/sd/a/ajross/holi-webjax/altzbin_mock$i/ --combine --nran 18 --boxsize 10000 --cellsize 10

done