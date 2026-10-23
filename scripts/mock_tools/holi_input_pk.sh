#!/bin/bash
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 4 clustering-stats --tracer LRG_input --zrange 0.5 0.85 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --stats_dir /pscratch/sd/a/ajross/holi-webjax/box10_mock$i/ --combine --nran 18 --boxsize 10000 --meshsize 700

 srun -N 1 -n 4 clustering-stats --tracer LRG_input --zrange 0.4 1.1 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --stats_dir /pscratch/sd/a/ajross/holi-webjax/box10_mock$i/ --combine --nran 18 --boxsize 10000 --meshsize 700
done