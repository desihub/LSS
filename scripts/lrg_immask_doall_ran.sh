#!/bin/bash

for i in `seq $1 $2`
do
    srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python readwrite_pixel_bitmask.py --tracer lrg --input $i --ran True
done

