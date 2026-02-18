#!/bin/bash
for ((i=$1;i<=$2;i++ ))
do
 srun -N 1 -n 4 clustering-stats --tracer LRG_input --zrange 0.5 0.85 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --stats_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --combine --nran 18 --analysis png_local

 srun -N 1 -n 4 clustering-stats --tracer LRG_input --zrange 0.4 1.1 --stats mesh2_spectrum --cat_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --stats_dir /pscratch/sd/a/ajross/holi-webjax/mock$i/ --combine --nran 18 --analysis png_local
done